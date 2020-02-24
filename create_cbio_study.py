'''
    for each section in yaml, perform:
        filtering ✔
        run cn data output generator ✔
        conversion to maf

    afterwards:
        merge mafs ✔
        merge gistic gene results ✔
        merge log seg data results ✔
        create cBio stidy meta files/connect to data files properly

    merge wgs with impact data?
    idea:   take two existing studies and compare files to see if any
            any overlap, take new description, name etc. as user
            input and merge exisitng non-overlapping data.
    '''

import click
import gzip
import yaml

from convert_vcf_to_maf import convert as convert_vcf_to_maf
from generate_outputs import extract, transform, load
from merge_outputs import merge_all_data as merge_outputs
from pathlib import Path


def merge_studies(path_to_output_study, path_to_external_study, output_dir):
    pass


def create_study(patient_yaml, path_to_output_study):
    '''
    data_CNA.txt, data_cna_hg19.seg, data_mutations.maf exist
    (in path_to_output_study)

    Please see cbioportal/docs/File-Formats.md on GitHub for examples
    '''

    meta_study = open(path_to_output_study + 'meta_study.txt', 'w+')
    
    type_of_cancer = click.prompt('Please enter type of cancer', default='ovary')
    cancer_study_identifier = click.prompt('Please enter a cancer study identifier', default='twins_shahlab_2020')
    name = click.prompt('Please enter a name', default='TWINS (Shah Lab, 2020)')
    description = click.prompt('Please enter a description', default='Mutation data in TWINS cases')
    short_name = click.prompt('Please enter a short name', default='TWINS (Shahlab)')
    
    meta_study.write('type_of_cancer: ' + type_of_cancer + '\n' \
                    + 'cancer_study_identifier: ' + cancer_study_identifier + '\n' \
                    + 'name: ' + name + '\n' \
                    + 'description: ' + description + '\n' \
                    + 'short_name: ' + short_name + '\n' \
                    + 'add_global_case_list: true\n')

    meta_mutations_extended = open(path_to_output_study + 'meta_mutations_extended.txt', 'w+')
    # TODO

    meta_CNA = open(path_to_output_study + 'meta_CNA.txt', 'w+')
    # TODO
    
    meta_cna_seg = open(path_to_output_study + 'meta_cna_seg.txt', 'w+')
    # TODO
    
    meta_clinical_sample = open(path_to_output_study + 'meta_clinical_sample.txt', 'w+')
    # TODO
    
    data_clinical_sample = open(path_to_output_study + 'data_clinical_sample.txt', 'w+')
    data_clinical_sample.write('#Patient Identifier\tSample Identifier\tCancer Type\n' \
                    + '#Patient Identifier\tSample Identifier\tCancer Type description\n' \
                    + '#STRING\tSTRING\tSTRING\n' \
                    + '#1\t1\t1\n' \
                    + '#PATIENT_ID\tSAMPLE_ID\tCANCER_TYPE\n')
    
    case_list_ids = []
    for patient, doc in patient_yaml.items():
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
                    + 'stable_id: ' + cancer_study_identifier + '_cna\n' \
                    + 'case_list_name: Sequenced samples profiled for mutations\n' \
                    + 'case_list_description: Sequenced samples profiled for mutations\n' \
                    + 'case_list_ids: ')

    cases_sequenced = open(path_to_output_study + 'case_lists/cases_sequenced.txt', 'w+')
    cases_sequenced.write('cancer_study_identifier: ' + cancer_study_identifier + '\n' \
                    + 'stable_id: ' + cancer_study_identifier + '_cna\n' \
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


def filter(sample_id, museq_vcf, strelka_vcf, work_dir):
    '''
    original code by Diljot Grewal

    required mofifications (done): take museq and strelka directly as input, output to temp_dir
    strip sample id from inputs and append to output filtered file
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
@click.option('--temp_dir', default='')
@click.option('--path_to_external_study')
@click.option('--path_to_merged_study')
def main(input_yaml, path_to_output_study, temp_dir, path_to_external_study):
    with open(input_yaml) as file:
        yaml_file = yaml.full_load(file)
        hgnc_file = yaml_file['id_mapping']
        gtf_file = yaml_file['gtf']
        
        for _, doc in yaml_file['patients'].items():
            for sample, doc in doc.items():
                museq_filtered = filter(sample, doc['museq_vcf'], doc['strelka_vcf'], temp_dir)
                
                with gzip.open(museq_filtered, 'rt') as museq_vcf:
                    # TODO: write converter
                    convert_vcf_to_maf(museq_vcf, sample, hgnc_file, temp_dir)
                
                with gzip.open(doc['titan_igv'], 'rt') as titan_igv, gzip.open(doc['titan_segs'], 'rt') as titan_segs:
                    generate_outputs(gtf_file, hgnc_file, titan_igv, titan_segs, sample, temp_dir)
                
        Path(path_to_output_study).mkdir(parents=True, exist_ok=True)
        merge_outputs(temp_dir, path_to_output_study)
        create_study(yaml_file['patients'], path_to_output_study)

    if path_to_external_study:
        merge_studies(path_to_output_study, path_to_external_study, path_to_merged_study)


if __name__ == '__main__':
    main()
