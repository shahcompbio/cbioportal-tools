'''
for each section in yaml, perform:
    filtering ✔
    run cn data output generator ✔
    conversion to maf ✔

afterwards:
    merge mafs ✔
    merge gistic gene results ✔
    merge log seg data results ✔
    create cBio stidy meta files/connect to data files properly ✔

'''

import click
import hmmcopy
import gzip
import yaml

from convert_vcf_to_maf import convert as convert_vcf_to_maf
from generate_outputs import extract, transform, load
from merge_outputs import merge_all_data as merge_outputs, merge_maf_data
from pathlib import Path


def create_study(study_info, path_to_output_study):
    '''
    prerequisites:
    data_CNA.txt, data_cna_hg19.seg, data_mutations_extended.maf exist
    (in path_to_output_study)

    Please see cbioportal/docs/File-Formats.md on GitHub for examples
    '''
    
    type_of_cancer = study_info['type_of_cancer']
    cancer_study_identifier = study_info['cancer_study_identifier']
    name = study_info['name']
    description = study_info['description']
    short_name = study_info['short_name']
    
    meta_study = open(path_to_output_study + 'meta_study.txt', 'w+')
    meta_study.write('type_of_cancer: ' + type_of_cancer + '\n' \
                    + 'cancer_study_identifier: ' + cancer_study_identifier + '\n' \
                    + 'name: ' + name + '\n' \
                    + 'description: ' + description + '\n' \
                    + 'short_name: ' + short_name + '\n' \
                    + 'add_global_case_list: true\n')

    meta_mutations_extended = open(path_to_output_study + 'meta_mutations_extended.txt', 'w+')
    meta_mutations_extended.write('cancer_study_identifier: ' + cancer_study_identifier + '\n' \
                    + 'genetic_alteration_type: MUTATION_EXTENDED\n' \
                    + 'datatype: MAF\n' \
                    + 'stable_id: mutations\n' \
                    + 'show_profile_in_analysis_tab: true\n' \
                    + 'profile_name: Mutations\n' \
                    + 'profile_description: Mutation data.\n' \
                    + 'data_filename: data_mutations_extended.maf\n')

    meta_CNA = open(path_to_output_study + 'meta_CNA.txt', 'w+')
    meta_CNA.write('cancer_study_identifier: ' + cancer_study_identifier + '\n' \
                    + 'genetic_alteration_type: COPY_NUMBER_ALTERATION\n' \
                    + 'datatype: DISCRETE\n' \
                    + 'stable_id: gistic\n' \
                    + 'show_profile_in_analysis_tab: true\n' \
                    + 'profile_name: Copy-number values\n' \
                    + 'profile_description: Copy-number values for each gene.\n' \
                    + 'data_filename: data_CNA.txt\n')
    
    meta_cna_seg = open(path_to_output_study + 'meta_cna_seg.txt', 'w+')
    meta_cna_seg.write('cancer_study_identifier: ' + cancer_study_identifier + '\n' \
                    + 'genetic_alteration_type: COPY_NUMBER_ALTERATION\n' \
                    + 'datatype: SEG\n' \
                    + 'reference_genome_id: hg19\n' \
                    + 'description: CNA data.\n' \
                    + 'data_filename: data_cna_hg19.seg\n')
    
    meta_clinical_sample = open(path_to_output_study + 'meta_clinical_sample.txt', 'w+')
    meta_clinical_sample.write('cancer_study_identifier: ' + cancer_study_identifier + '\n' \
                    + 'genetic_alteration_type: CLINICAL\n' \
                    + 'datatype: SAMPLE_ATTRIBUTES\n' \
                    + 'data_filename: data_clinical_sample.txt\n')
    
    data_clinical_sample = open(path_to_output_study + 'data_clinical_sample.txt', 'w+')
    data_clinical_sample.write('#Patient Identifier\tSample Identifier\tCancer Type\n' \
                    + '#Patient Identifier\tSample Identifier\tCancer Type description\n' \
                    + '#STRING\tSTRING\tSTRING\n' \
                    + '#1\t1\t1\n' \
                    + 'PATIENT_ID\tSAMPLE_ID\tCANCER_TYPE\n')
    
    case_list_ids = []
    for patient, doc in study_info['patients'].items():
        for sample, _ in doc.items():
            data_clinical_sample.write(patient + '\t' + sample + '\t' + type_of_cancer.upper() + '\n')
            case_list_ids.append(sample)
    
    Path(path_to_output_study + 'case_lists/').mkdir(parents=True, exist_ok=True)

    cases_cna = open(path_to_output_study + 'case_lists/cases_cna.txt', 'w+')
    cases_cna.write('cancer_study_identifier: ' + cancer_study_identifier + '\n' \
                    + 'stable_id: ' + cancer_study_identifier + '_cna\n' \
                    + 'case_list_name: Samples profiled for mutations\n' \
                    + 'case_list_description: Samples profiled for mutations\n' \
                    + 'case_list_ids: ')

    cases_cnaseq = open(path_to_output_study + 'case_lists/cases_cnaseq.txt', 'w+')
    cases_cnaseq.write('cancer_study_identifier: ' + cancer_study_identifier + '\n' \
                    + 'stable_id: ' + cancer_study_identifier + '_cnaseq\n' \
                    + 'case_list_name: Sequenced samples profiled for mutations\n' \
                    + 'case_list_description: Sequenced samples profiled for mutations\n' \
                    + 'case_list_ids: ')

    cases_sequenced = open(path_to_output_study + 'case_lists/cases_sequenced.txt', 'w+')
    cases_sequenced.write('cancer_study_identifier: ' + cancer_study_identifier + '\n' \
                    + 'stable_id: ' + cancer_study_identifier + '_sequenced\n' \
                    + 'case_list_name: All sequenced samples\n' \
                    + 'case_list_description: All sequenced samples\n' \
                    + 'case_list_ids: ')

    for sample_id in case_list_ids[:-1]:
        cases_cna.write(sample_id + '\t')
        cases_cnaseq.write(sample_id + '\t')
        cases_sequenced.write(sample_id + '\t')
    
    for file in [cases_cna, cases_cnaseq, cases_sequenced]:
        file.write(case_list_ids[-1] + '\n')


def generate_outputs(gtf_file, hgnc_file, titan_igv, titan_segs, sample_id, output_dir): 
    extracted_file = extract(gtf_file, hgnc_file, titan_igv, titan_segs)
    gene_dict, seg_dict = transform(extracted_file, show_missing_hugo=False, show_missing_entrez=False, show_missing_both=False)
    load(gene_dict, seg_dict, sample_id, output_dir, output_gistic_gene=True, output_integer_gene=False, output_log_seg=True, output_integer_seg=False)


def filter_vcfs(sample_id, museq_vcf, strelka_vcf, work_dir):
    '''
    original code by Diljot Grewal

    museq_paired and strekla_snv,
    take position intersection plus probability filter of 0.85
    (keep positions >= 0.85 that are in both)

    modifications: take museq and strelka directly as input, output
    to temp_dir, take sample id as input and append to output
    filtered filename
    '''

    museq_filtered = work_dir + '{}_museq_filtered.vcf.gz'.format(sample_id)

    strelka_ref = set()

    with gzip.open(strelka_vcf, 'rt') as strelka_data:
        for line in strelka_data:
            if line.startswith('#'):
                continue

            line = line.strip().split()
            
            chrom = line[0]
            pos = line[1]
            
            strelka_ref.add((chrom, pos))

    with gzip.open(museq_vcf, 'rt') as museq_data, gzip.open(museq_filtered, 'wt') as museqout:
        for line in museq_data:
            if line.startswith('#'):
                museqout.write(line)
                continue
            
            line = line.strip().split()
            chrom = line[0]
            pos = line[1]
            
            if ((chrom, pos))  not in strelka_ref:
                continue
            
            pr = line[7].split(';')[0].split('=')[1]
            if float(pr) < 0.85:
                continue
            
            outstr = '\t'.join(line)+'\n'
            museqout.write(outstr)

    return museq_filtered


@click.command()
@click.argument('input_yaml')
@click.argument('path_to_output_study')
@click.argument('temp_dir')
def main(input_yaml, path_to_output_study, temp_dir):
    if not path_to_output_study.endswith('/'):
        path_to_output_study = path_to_output_study + '/'
    if not temp_dir.endswith('/'):
        temp_dir = temp_dir + '/'

    Path(path_to_output_study).mkdir(parents=True, exist_ok=True)
    Path(temp_dir).mkdir(parents=True, exist_ok=True)

    with open(input_yaml) as file:
        yaml_file = yaml.full_load(file)
        hgnc_file = yaml_file['id_mapping']
        gtf_file = yaml_file['gtf']
        
        create_study(yaml_file, path_to_output_study)

        for patient_id, doc in yaml_file['patients'].items():
            for sample, doc in doc.items():
                museq_filtered = filter_vcfs(sample, doc['museq_vcf'], doc['strelka_vcf'], temp_dir)

                dataset_id = f'{patient_id}-{sample}-snvs'
                convert_vcf_to_maf(museq_filtered, sample, temp_dir, dataset_id)

                dataset_id = f'{patient_id}-{sample}-indels'
                convert_vcf_to_maf(doc['strelka_indel_vcf'], sample, temp_dir, dataset_id)

                if doc['datatype'] == 'WGS':
                    with gzip.open(doc['titan_segs'], 'rt') as titan_segs:
                        generate_outputs(gtf_file, hgnc_file, doc['titan_igv'], titan_segs, sample, temp_dir)        

                elif doc['datatype'] == 'SCWGS':
                    cnv = hmmcopy.read_copy_data(doc['hmmcopy_csv'], filter_normal=doc['filter_normal'])
                    genes = hmmcopy.read_gene_data(gtf_file)
                    
                    overlapping = hmmcopy.calculate_gene_copy(cnv, genes)
                    hmmcopy.convert_to_transform_format(overlapping, hgnc_file, temp_dir)

                    hmmcopy_extract = open(temp_dir + 'hmmcopy_extract', 'r')
                    gene_dict, seg_dict = transform(hmmcopy_extract, show_missing_hugo=False, show_missing_entrez=False, show_missing_both=False)
                    load(gene_dict, seg_dict, sample, temp_dir, output_gistic_gene=True, output_integer_gene=False, output_log_seg=True, output_integer_seg=False)

                else:
                    raise ValueError(f'unrecognized data type {doc["datatype"]}')

    merge_outputs(temp_dir, path_to_output_study)


if __name__ == '__main__':
    main()
