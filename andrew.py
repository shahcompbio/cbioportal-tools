import click
import pandas as pd
import yaml

from create_cbio_study import create_study
from pathlib import Path


def load_remixt(sample, filename, min_ploidy=None, max_ploidy=None):
    cn_data = {}
    stats_data = []
    
    with pd.HDFStore(filename) as store:
        stats = store['stats']
        stats = stats[stats['proportion_divergent'] < 0.5]
        
        if max_ploidy:
            stats = stats[stats['ploidy'] < max_ploidy]
        
        if min_ploidy:
            stats = stats[stats['ploidy'] > min_ploidy]
        
        stats = stats.sort_values('elbo').iloc[-1]
        stats['sample'] = sample

        init_id = stats['init_id']

        cn = store[f'/solutions/solution_{init_id}/cn']
        cn['segment_length'] = cn['end'] - cn['start'] + 1
        cn['length_ratio'] = cn['length'] / cn['segment_length']

        mix = store[f'/solutions/solution_{init_id}/mix']

        stats['normal_proportion'] = mix[0]

        cn_data[sample] = cn
        stats_data.append(stats)

    stats_data = pd.DataFrame(stats_data)
    stats_data['tumour_proportion'] = 1. - stats_data['normal_proportion']

    stats_data[[
        'sample', 'ploidy', 'proportion_divergent',
        'tumour_proportion', 'proportion_divergent', 'elbo']].sort_values('sample').to_csv('remixt_stats.csv', index=False)

    stats_data[[
        'sample', 'ploidy', 'proportion_divergent',
        'tumour_proportion', 'proportion_divergent', 'elbo']].sort_values('sample')
    
    return cn_data, stats_data


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
    hgnc = pd.read_csv(hgnc_file, delimiter='\t', dtype=str)
    hgnc = hgnc.rename(columns={'Approved symbol': 'Hugo_Symbol', 'NCBI Gene ID': 'Entrez_Gene_Id', 'Ensembl gene ID': 'gene_id'})
    genes = genes.merge(hgnc_file, on=['gene_id'], how='left')

    return genes


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
        genes_with_ids = hgnc_lookup(genes, hgnc_file)
        # TODO

        for patient_id, patient_data in yaml_file['patients'].items():
            for sample, sample_data in patient_data.items():
                if sample_data['datatype'] == 'WGS':
                    pass
                    # TODO 

    merge_outputs(temp_dir, path_to_output_study)
