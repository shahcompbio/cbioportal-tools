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
import numpy as np
import pandas as pd
import requests
import wgs_analysis.algorithms.cnv
import yaml

from convert_vcf_to_maf import get_vcf_data, write_vcf, write_allele_counts, generate_mafs
from generate_outputs import transform, load
from merge_outputs import merge_all_data as merge_outputs
from os import path
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
    
    maf.to_csv(temp_dir + sample_id + '.maf', index=None, sep='\t')


def read_gene_data(gtf):
    data = pd.read_csv(
        gtf,
        delimiter='\t',
        names=['chromosome', 'gene_start', 'gene_end', 'info'],
        usecols=[0,3,4,8],
        converters={'chromosome': str},
    )

    def extract_info(info):
        info_dict = {}
        for a in info.split('; '):
            k, v = a.split(' ')
            info_dict[k] = v.strip(';').strip('"')
        return info_dict
    
    data['info'] = data['info'].apply(extract_info)
    data['gene_id'] = data['info'].apply(lambda a: a['gene_id'])
    data['gene_name'] = data['info'].apply(lambda a: a['gene_name'])

    data = data.groupby(['chromosome', 'gene_id', 'gene_name']).agg({'gene_start':'min', 'gene_end':'max'}).reset_index()

    return data


def hgnc_lookup(genes, hgnc_file):
    genes_page_0 = requests.get('https://www.cbioportal.org/api/genes')
    genes_page_1 = requests.get('https://www.cbioportal.org/api/genes?pageNumber=1')
    gene_request = genes_page_0.json() + genes_page_1.json()

    cbio_genes = pd.DataFrame(gene_request, dtype=str)
    cbio_genes.drop('type', axis=1, inplace=True)
    cbio_genes.rename(columns={'hugoGeneSymbol': 'Hugo_Symbol', 'entrezGeneId': 'Entrez_Gene_Id'}, inplace=True)
    cbio_genes['Hugo_Symbol'] = cbio_genes['Hugo_Symbol'].str.upper()

    hgnc = pd.read_csv(hgnc_file, delimiter='\t', dtype=str)
    
    hgnc.dropna(subset=['Ensembl gene ID', 'Ensembl ID(supplied by Ensembl)'], how='all', inplace=True)
    hgnc.loc[hgnc['Ensembl gene ID'].isna(), 'Ensembl gene ID'] = hgnc['Ensembl ID(supplied by Ensembl)']
    hgnc.drop('Ensembl ID(supplied by Ensembl)', axis=1, inplace=True)
    hgnc.rename(columns={'Approved symbol': 'Hugo_Symbol', 'Ensembl gene ID': 'gene_id'}, inplace=True)

    final_genes = genes.merge(hgnc, on=['gene_id'], how='left')
    final_genes['Hugo_Symbol'] = final_genes['Hugo_Symbol'].str.upper()
    final_genes.dropna(subset=['Hugo_Symbol'], inplace=True)

    # hugo_not_in_cbio
    hugo_not_in_cbio = hgnc.merge(cbio_genes, on=['Hugo_Symbol'], how='left')
    hugo_not_in_cbio = hugo_not_in_cbio[hugo_not_in_cbio['Entrez_Gene_Id'].isna()]
    hugo_not_in_cbio.drop('Entrez_Gene_Id', axis=1, inplace=True)
    hugo_not_in_cbio.to_csv('counts/hugo_not_in_cbio.txt', index=None, sep='\t')
    cbio_counts = hugo_not_in_cbio.groupby(['Locus group', 'Locus type']).size()
    cbio_counts = cbio_counts.reset_index()
    cbio_counts.rename(columns={0: 'Count'}, inplace=True)
    cbio_counts.to_csv('counts/hugo_not_in_cbio_counts.txt', index=None, sep='\t')

    # cbio_not_in_hugo
    cbio_not_in_hugo = hgnc.merge(cbio_genes, on=['Hugo_Symbol'], how='right')
    cbio_not_in_hugo = cbio_not_in_hugo[cbio_not_in_hugo['gene_id'].isna()]
    cbio_not_in_hugo = cbio_not_in_hugo[['Hugo_Symbol', 'Entrez_Gene_Id']]
    cbio_not_in_hugo.to_csv('counts/cbio_not_in_hugo.txt', index=None, sep='\t')

    # hugo_not_in_gtf
    hugo_not_in_gtf = final_genes[['sample', 'gene_id']].merge(hgnc, on=['gene_id'], how='right')
    hugo_not_in_gtf = hugo_not_in_gtf[hugo_not_in_gtf['sample'].isna()]
    hugo_not_in_gtf.drop(['sample'], axis=1, inplace=True)
    hugo_not_in_gtf.to_csv('counts/hugo_not_in_gtf.txt', index=None, sep='\t')
    gtf_counts = hugo_not_in_gtf.groupby(['Locus group', 'Locus type']).size()
    gtf_counts = gtf_counts.reset_index()
    gtf_counts.rename(columns={0: 'Count'}, inplace=True)
    gtf_counts.to_csv('counts/hugo_not_in_gtf_counts.txt', index=None, sep='\t')
    
    final_genes = final_genes.merge(cbio_genes, on=['Hugo_Symbol'], how='left')
    
    # manual additions as per request
    final_genes.loc[final_genes['gene_id'] == 'ENSG00000133706', 'Hugo_Symbol'] = 'LARS'
    final_genes.loc[final_genes['gene_id'] == 'ENSG00000237452', 'Hugo_Symbol'] = 'BHMG1'

    return final_genes


def process_remixt(sample, sample_data):
    with pd.HDFStore(sample_data['remixt']) as store:
        stats = store['stats']
        stats = stats[stats['proportion_divergent'] < 0.5]
        
        if 'max_ploidy' in sample_data:
            stats = stats[stats['ploidy'] < max_ploidy]
        
        if 'min_ploidy' in sample_data:
            stats = stats[stats['ploidy'] > min_ploidy]
        
        stats = stats.sort_values('elbo').iloc[-1]
        stats['sample'] = sample

        init_id = stats['init_id']

        cn = store[f'/solutions/solution_{init_id}/cn']
        cn['segment_length'] = cn['end'] - cn['start'] + 1
        cn['length_ratio'] = cn['length'] / cn['segment_length']

        mix = store[f'/solutions/solution_{init_id}/mix']

        stats['normal_proportion'] = mix[0]

        return cn, stats


def generate_aggregated_cn(cn_data):
    aggregated_cn_data = {}

    for sample, cn in cn_data.items():
        aggregated_cn_data[sample] = wgs_analysis.algorithms.cnv.aggregate_adjacent(
            cn,
            value_cols=['major_0', 'minor_0', 'major_1', 'minor_1', 'major_2', 'minor_2'],
            stable_cols=['major_0', 'minor_0', 'major_1', 'minor_1', 'major_2', 'minor_2'],
            length_normalized_cols=['major_raw', 'minor_raw'],
        )

    return aggregated_cn_data


def generate_genes_cn(aggregated_cn_data, genes):
    genes_cn_data = {}

    for sample, cn in aggregated_cn_data.items():
        cn['width'] = cn['end'] - cn['start']
        cn['total_raw'] = cn['major_raw'] + cn['minor_raw']

        genes_cn_data[sample] = wgs_analysis.algorithms.cnv.calculate_gene_copy(
            cn, genes,
            [
                'major_raw',
                'minor_raw',
                'total_raw',
                'major_1',
                'minor_1',
                'major_2',
                'minor_2',
            ])

    return genes_cn_data


def generate_amp(genes_cn_data, stats_data, genes):
    amp_data = []

    for sample, data in genes_cn_data.items():
        normalize = (
            data
            .groupby('gene_id')['overlap_width']
            .sum().rename('sum_overlap_width').reset_index())

        data['total_raw_weighted'] = data['total_raw'] * data['overlap_width']

        data = data.groupby(['gene_id'])['total_raw_weighted'].sum().reset_index()
        data = data.merge(normalize)
        data['total_raw_mean'] = data['total_raw_weighted'] / data['sum_overlap_width']
        data['sample'] = sample

        amp_data.append(data[['gene_id', 'total_raw_mean', 'sample']])

    amp_data = pd.concat(amp_data)

    gene_cols = [
        'gene_id',
        'chromosome',
        'gene_start',
        'gene_end',
        'gene_name',
    ]
    amp_data = amp_data.merge(genes[gene_cols])
    amp_data = amp_data.merge(stats_data[['sample', 'ploidy']])
    amp_data['log_change'] = np.log2(amp_data['total_raw_mean'] / amp_data['ploidy'])
    amp_data['log_change'] = amp_data['log_change'].fillna(np.exp(-8))
    amp_data.loc[amp_data['log_change'] == np.NINF, 'log_change'] = np.exp(-8)

    amp_data = amp_data[['gene_id', 'gene_name', 'sample', 'log_change']]

    return amp_data


def generate_hdel(genes_cn_data, genes):
    hdel_data = []

    for sample, data in genes_cn_data.items():
        data = data[data['total_raw'] < 0.5]
        data = data.groupby(['gene_id'])['overlap_width'].sum().rename('hdel_width').reset_index()
        data = data[data['hdel_width'] > 10000]
        data['sample'] = sample

        hdel_data.append(data[['gene_id', 'hdel_width', 'sample']])

    hdel_data = pd.concat(hdel_data)

    gene_cols = [
        'gene_id',
        'chromosome',
        'gene_start',
        'gene_end',
        'gene_name',
    ]
    hdel_data = hdel_data.merge(genes[gene_cols])

    hdel_data = hdel_data[['gene_id', 'sample']]

    return hdel_data


def generate_gistic_outputs(gistic_data, hdel_data, path_to_output_study, hgnc_file):
    # Classify by log change
    gistic_data['gistic_value'] = 2
    gistic_data.loc[gistic_data['log_change'] < 1, 'gistic_value'] = 1
    gistic_data.loc[gistic_data['log_change'] < 0.5, 'gistic_value'] = 0
    gistic_data.loc[gistic_data['log_change'] < -0.5, 'gistic_value'] = -1
    
    # Merge hdels
    hdel_data['is_hdel'] = 1
    gistic_data = gistic_data.merge(hdel_data[['gene_id', 'sample', 'is_hdel']], how='left')
    gistic_data['is_hdel'] = gistic_data['is_hdel'].fillna(0).astype(int)
    gistic_data.loc[gistic_data['is_hdel'] == 1, 'gistic_value'] = -2

    # Testing gistic_data generation
    gistic_data = hgnc_lookup(gistic_data, hgnc_file)
    gistic_data = gistic_data[['Hugo_Symbol', 'Entrez_Gene_Id', 'sample', 'gistic_value']]
    gistic_matrix = gistic_data.set_index(['Hugo_Symbol', 'Entrez_Gene_Id', 'sample'])['gistic_value'].unstack()
    gistic_matrix.reset_index(inplace=True)
    gistic_matrix.to_csv(path_to_output_study + 'data_CNA.txt', index=None, sep='\t')


def generate_seg_outputs(aggregated_cn_data, temp_dir):
    # clean up segs and write to disk
    for sample in aggregated_cn_data:
        aggregated_cn_data[sample]['sample'] = sample
        aggregated_cn_data[sample] = aggregated_cn_data[sample].merge(stats_data[['sample', 'ploidy']])
        aggregated_cn_data[sample]['total_raw'] = aggregated_cn_data[sample]['major_raw'] + aggregated_cn_data[sample]['minor_raw']
        aggregated_cn_data[sample]['seg.mean'] = np.log2(aggregated_cn_data[sample]['total_raw'] / aggregated_cn_data[sample]['ploidy'])
        aggregated_cn_data[sample]['num.mark'] = (aggregated_cn_data[sample]['length'] / 500000).astype(int)
        aggregated_cn_data[sample] = aggregated_cn_data[sample].rename(columns={'sample': 'ID', 'chromosome': 'chrom', 'start': 'loc.start', 'end': 'loc.end'})
        aggregated_cn_data[sample] = aggregated_cn_data[sample][['ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean']]
        aggregated_cn_data[sample]['seg.mean'] = aggregated_cn_data[sample]['seg.mean'].fillna(np.exp(-8))
        aggregated_cn_data[sample].loc[aggregated_cn_data[sample]['seg.mean'] == np.NINF, 'seg.mean'] = np.exp(-8)
        aggregated_cn_data[sample].to_csv(temp_dir + sample + '_log_seg_data.seg', index=None, sep='\t')


def clean_up_stats(stats_data):
    stats_data = pd.DataFrame(stats_data)
    stats_data['tumour_proportion'] = 1. - stats_data['normal_proportion']

    stats_data = stats_data[[
        'sample', 'ploidy', 'proportion_divergent',
        'tumour_proportion', 'proportion_divergent', 'elbo']].sort_values('sample')

    return stats_data

        
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
        genes = read_gene_data(gtf_file)

        cn_data = {}
        stats_data = []
        for patient_id, patient_data in yaml_file['patients'].items():
            for sample, sample_data in patient_data.items():
                if sample_data['datatype'] == 'WGS':
                    vcf_files[sample] = []

                    if 'museq_vcf' in sample_data and 'strelka_vcf' in sample_data:
                        museq_filtered = filter_vcfs(sample, sample_data['museq_vcf'], sample_data['strelka_vcf'], temp_dir)
                        vcf_files[sample].append(museq_filtered)
                    
                    if 'strelka_indel_vcf' in sample_data:
                        vcf_files[sample].append(sample_data['strelka_indel_vcf'])

                    if 'remixt' in sample_data:
                        cn, stats = process_remixt(sample, sample_data)
                        cn_data[sample] = cn
                        stats_data.append(stats)            

                elif sample_data['datatype'] == 'SCWGS':
                    hmmcopy_list = []
                    snv_counts = []
                    vcf_files[sample] = []

                    for library_id, library_data in sample_data.items():
                        if library_id == 'datatype':
                            continue

                        if 'museq_vcf' in library_data:
                            vcf_files[sample].append(library_data['museq_vcf'])
                        
                        if 'strelka_vcf' in library_data:
                            vcf_files[sample].append(library_data['strelka_vcf'])
                        
                        if 'strelka_indel_vcf' in library_data:
                            vcf_files[sample].append(library_data['strelka_indel_vcf'])

                        hmmcopy_list.append(library_data['hmmcopy_csv'])  
                        
                        if 'snv_counts_csv' in library_data:
                            snv_counts.append(library_data['snv_counts_csv'])

                    merge_hmmcopy(hmmcopy_list, temp_dir)
                    
                    if snv_counts:
                        calculate_counts(snv_counts, patient_id, sample, temp_dir)
                        t_file = Path(temp_dir + sample + '_tumour_counts.csv')
                    
                    cnv = hmmcopy.read_copy_data(temp_dir + 'hmmcopy_csv', filter_normal=False)
                    hmmcopy_genes = hmmcopy.read_gene_data(gtf_file)
                    
                    overlapping = hmmcopy.calculate_gene_copy(cnv, hmmcopy_genes)
                    hmmcopy.convert_to_transform_format(overlapping, hgnc_file, temp_dir)

                    hmmcopy_extract = open(temp_dir + 'hmmcopy_extract', 'r')
                    gene_dict, seg_dict = transform(hmmcopy_extract, show_missing_hugo=False, show_missing_entrez=False, show_missing_both=False)
                    load(gene_dict, seg_dict, sample, temp_dir, output_gistic_gene=True, output_integer_gene=False, output_log_seg=True, output_integer_seg=False)

                else:
                    raise ValueError(f'unrecognized data type {sample_data["datatype"]}')
        
        stats_data = clean_up_stats(stats_data)

        aggregated_cn_data = generate_aggregated_cn(cn_data)
        genes_cn_data = generate_genes_cn(aggregated_cn_data, genes)
        amp_data = generate_amp(genes_cn_data, stats_data, genes)
        hdel_data = generate_hdel(genes_cn_data, genes)

        generate_gistic_outputs(amp_data, hdel_data, path_to_output_study, hgnc_file)
        generate_seg_outputs(aggregated_cn_data, temp_dir)

        vcf_outputs = {sample: path.join(temp_dir, '{}.vcf'.format(sample)) for sample in vcf_files}
        csv_outputs = {sample: path.join(temp_dir, '{}.csv'.format(sample)) for sample in vcf_files}
        maf_outputs = {sample: path.join(temp_dir, '{}.maf'.format(sample)) for sample in vcf_files}

        for sample, sample_vcf_files in vcf_files.items():
            vcfdata = get_vcf_data(sample_vcf_files)
            write_vcf(vcf_outputs[sample], sample_vcf_files, vcfdata, temp_dir, sample)
            write_allele_counts(csv_outputs[sample], vcfdata)

        generate_mafs(vcf_files, temp_dir, maf_outputs, vcf_outputs)

        for patient_id, patient_data in yaml_file['patients'].items():
            for sample, _ in patient_data.items():
                n_file = Path(temp_dir + sample + '.csv')
            
                maf = Path(temp_dir + sample + '.maf')
                if n_file.is_file() and maf.is_file():
                    add_counts_to_maf(patient_id, sample, temp_dir)

    merge_outputs(temp_dir, path_to_output_study)


if __name__ == '__main__':
    main()
