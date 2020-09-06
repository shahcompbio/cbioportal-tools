import numpy as np
import pandas as pd
import wgs_analysis.algorithms.cnv

from utils import hgnc_lookup


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


def load_data(sample, sample_data):
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
    aggregated_cn_data = wgs_analysis.algorithms.cnv.aggregate_adjacent(
        cn,
        value_cols=['major_0', 'minor_0', 'major_1', 'minor_1', 'major_2', 'minor_2'],
        stable_cols=['major_0', 'minor_0', 'major_1', 'minor_1', 'major_2', 'minor_2'],
        length_normalized_cols=['major_raw', 'minor_raw'],
    )

    return aggregated_cn_data


def generate_genes_cn(aggregated_cn_data, genes):
    genes_cn_data = wgs_analysis.algorithms.cnv.calculate_gene_copy(
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
    amp_data = genes_cn_data.copy()

    normalize = (
        amp_data
        .groupby('gene_id')['overlap_width']
        .sum().rename('sum_overlap_width').reset_index())

    amp_data['total_raw_weighted'] = amp_data['total_raw'] * amp_data['overlap_width']

    amp_data = amp_data.groupby(['gene_id'])['total_raw_weighted'].sum().reset_index()
    amp_data = amp_data.merge(normalize)
    amp_data['total_raw_mean'] = amp_data['total_raw_weighted'] / amp_data['sum_overlap_width']

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

    amp_data = amp_data[['gene_id', 'gene_name', 'log_change']]

    return amp_data


def generate_hdel(genes_cn_data, genes):
    hdel_data = genes_cn_data.copy()

    hdel_data = hdel_data[hdel_data['total_raw'] < 0.5]
    hdel_data = hdel_data.groupby(['gene_id'])['overlap_width'].sum().rename('hdel_width').reset_index()
    hdel_data = hdel_data[hdel_data['hdel_width'] > 10000]

    gene_cols = [
        'gene_id',
        'chromosome',
        'gene_start',
        'gene_end',
        'gene_name',
    ]
    hdel_data = hdel_data.merge(genes[gene_cols])

    hdel_data = hdel_data[['gene_id']]

    return hdel_data


def generate_gistic_outputs(filename, gistic_data, hdel_data, hgnc_file):
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

    # Gistic_data generation
    gistic_data = hgnc_lookup(gistic_data, hgnc_file)
    gistic_data = gistic_data[['Hugo_Symbol', 'Entrez_Gene_Id', 'sample', 'gistic_value']]
    gistic_matrix = gistic_data.set_index(['Hugo_Symbol', 'Entrez_Gene_Id', 'sample'])['gistic_value'].unstack()
    gistic_matrix.reset_index(inplace=True)
    gistic_matrix.to_csv(filename, index=None, sep='\t')


def generate_seg_outputs(filename, aggregated_cn_data, stats_data):
    # Clean up segs and write to disk
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
        aggregated_cn_data[sample].to_csv(filename, index=None, sep='\t')


def clean_up_stats(stats_data):
    stats_data = pd.DataFrame(stats_data)
    stats_data['tumour_proportion'] = 1. - stats_data['normal_proportion']

    stats_data = stats_data[[
        'ploidy', 'proportion_divergent',
        'tumour_proportion', 'proportion_divergent', 'elbo']].sort_values('sample')

    return stats_data
