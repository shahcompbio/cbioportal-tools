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
import pandas as pd
import subprocess
import yaml

from convert_vcf_to_maf import convert as convert_vcf_to_maf
from generate_outputs import extract, transform, load
from merge_outputs import merge_all_data as merge_outputs
from os import path
from pathlib import Path
from run import get_vcf_data, write_vcf, write_allele_counts, generate_mafs


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
    for patient, patient_data in study_info['patients'].items():
        for sample, _ in patient_data.items():
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


def merge_hmmcopy(hmmcopy_files, temp_dir):
    final_df = pd.read_csv(hmmcopy_files.pop(0), dtype=object)
    
    for file in hmmcopy_files:
        df = pd.read_csv(file, dtype=object)
        final_df = pd.concat([df, final_df], axis=0, ignore_index=True)
    
    final_df.to_csv(temp_dir + 'hmmcopy_csv', index=None)


def calculate_counts(counts_files, patient_id, sample_id, temp_dir):
    usecols = ['chrom','coord','ref','alt', 'ref_counts', 'alt_counts']
    final_df = pd.DataFrame(columns=usecols)
    for counts_file in counts_files:
        for df in pd.read_csv(counts_file, chunksize=1e6, usecols=usecols):
            final_df = pd.concat([df, final_df], axis=0, ignore_index=True)
            final_df = final_df.groupby(['chrom','coord','ref','alt'], as_index=False).agg('sum')

    final_df = final_df.rename(columns={'ref_counts': 't_ref_count', 'alt_counts': 't_alt_count'})
    final_df.to_csv(temp_dir + sample_id + '_tumour_counts.csv', index=None, sep='\t')


def add_counts_to_maf(patient_id, sample_id, temp_dir):
    n_counts = pd.read_csv(temp_dir + sample_id + '.csv', dtype=object, sep='\t').rename(columns={'chrom': 'Chromosome', 'coord': 'Start_Position', 'ref': 'Reference_Allele', 'alt': 'Tumor_Seq_Allele2'})
    
    maf = pd.read_csv(temp_dir + sample_id + '.maf', dtype=object, sep='\t', skiprows=1)
    maf = maf.drop(columns=['n_ref_count', 'n_alt_count'])
    maf = maf.merge(n_counts, on=['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2'], how='left')
    
    if path.exists(temp_dir + sample_id + '_tumour_counts.csv'):
        t_counts = pd.read_csv(temp_dir + sample_id + '_tumour_counts.csv', dtype=object, sep='\t').rename(columns={'chrom': 'Chromosome', 'coord': 'Start_Position', 'ref': 'Reference_Allele', 'alt': 'Tumor_Seq_Allele2'})
        maf = maf.drop(columns=['t_ref_count', 't_alt_count'])
        maf = maf.merge(t_counts, on=['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2'], how='left')
    
    maf.to_csv(temp_dir + sample_id + '-generated.maf', index=None, sep='\t')


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
        vcf_files = {}
        
        create_study(yaml_file, path_to_output_study)

        for patient_id, patient_data in yaml_file['patients'].items():
            for sample, sample_data in patient_data.items():
                if sample_data['datatype'] == 'WGS':
                    if 'museq_vcf' in sample_data and 'strelka_vcf' in sample_data:
                        museq_filtered = filter_vcfs(sample, sample_data['museq_vcf'], sample_data['strelka_vcf'], temp_dir)
                        dataset_id = f'{patient_id}-{sample}-snvs'
                        convert_vcf_to_maf(museq_filtered, sample, dataset_id, temp_dir)
                    
                    if 'strelka_indel_vcf' in sample_data:
                        dataset_id = f'{patient_id}-{sample}-indels'
                        convert_vcf_to_maf(sample_data['strelka_indel_vcf'], sample, dataset_id, temp_dir)

                    with gzip.open(sample_data['titan_segs'], 'rt') as titan_segs:
                        generate_outputs(gtf_file, hgnc_file, sample_data['titan_igv'], titan_segs, sample, temp_dir)        

                elif sample_data['datatype'] == 'SCWGS':
                    hmmcopy_list = []
                    snv_counts = []
                    vcf_files[sample] = []

                    for library_id, library_data in sample_data.items():
                        if library_id == 'datatype':
                            continue

                        vcf_files[sample].append(library_data['museq_vcf'])
                        vcf_files[sample].append(library_data['strelka_vcf'])
                        vcf_files[sample].append(library_data['strelka_indel_vcf'])

                        hmmcopy_list.append(sample_data['hmmcopy_csv'])  
                        
                        if 'snv_counts_csv' in sample_data:
                            snv_counts.append(sample_data['snv_counts_csv'])

                    merge_hmmcopy(hmmcopy_list, temp_dir)
                    
                    if snv_counts:
                        calculate_counts(snv_counts, patient_id, sample, temp_dir)
                        t_file = Path(temp_dir + sample + '_tumour_counts.csv')
                    
                    cnv = hmmcopy.read_copy_data(temp_dir + 'hmmcopy_csv', filter_normal=False)
                    genes = hmmcopy.read_gene_data(gtf_file)
                    
                    overlapping = hmmcopy.calculate_gene_copy(cnv, genes)
                    hmmcopy.convert_to_transform_format(overlapping, hgnc_file, temp_dir)

                    hmmcopy_extract = open(temp_dir + 'hmmcopy_extract', 'r')
                    gene_dict, seg_dict = transform(hmmcopy_extract, show_missing_hugo=False, show_missing_entrez=False, show_missing_both=False)
                    load(gene_dict, seg_dict, sample, temp_dir, output_gistic_gene=True, output_integer_gene=False, output_log_seg=True, output_integer_seg=False)

                else:
                    raise ValueError(f'unrecognized data type {sample_data["datatype"]}')

        vcf_outputs = {sample: os.path.join(temp_dir, '{}.vcf'.format(sample)) for sample in vcf_files}
        csv_ouputs = {sample: os.path.join(temp_dir, '{}.csv'.format(sample)) for sample in vcf_files}
        maf_outputs = {sample: os.path.join(temp_dir, '{}.maf'.format(sample)) for sample in vcf_files}

        vep_dir = '/juno/work/shah/svatrt/vcf2maf/cache/'

        for sample, sample_vcf_files in vcf_files.items():
            vcfdata = get_vcf_data(sample_vcf_files)
            write_vcf(vcf_outputs[sample], sample_vcf_files, vcfdata, temp_dir, sample)
            write_allele_counts(csv_outputs[sample], vcfdata)

        generate_mafs(vcf_files, vep_dir, temp_dir, maf_outputs, vcf_outputs)

        for sample in vcf_files:
            n_file = Path(temp_dir + sample + '.csv')
            
            maf = Path(temp_dir + sample + '.maf')
            if n_file.is_file() and maf.is_file():
                add_counts_to_maf(patient_id, sample, temp_dir)

    merge_outputs(temp_dir, path_to_output_study)


if __name__ == '__main__':
    main()
