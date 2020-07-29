import numpy as np
import pandas as pd
import requests
import wgs_analysis.algorithms.cnv


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

    # Generate hugo_not_in_cbio info
    hugo_not_in_cbio = hgnc.merge(cbio_genes, on=['Hugo_Symbol'], how='left')
    hugo_not_in_cbio = hugo_not_in_cbio[hugo_not_in_cbio['Entrez_Gene_Id'].isna()]
    hugo_not_in_cbio.drop('Entrez_Gene_Id', axis=1, inplace=True)
    hugo_not_in_cbio.to_csv('counts/hugo_not_in_cbio.txt', index=None, sep='\t')
    cbio_counts = hugo_not_in_cbio.groupby(['Locus group', 'Locus type']).size()
    cbio_counts = cbio_counts.reset_index()
    cbio_counts.rename(columns={0: 'Count'}, inplace=True)
    cbio_counts.to_csv('counts/hugo_not_in_cbio_counts.txt', index=None, sep='\t')

    # Generate cbio_not_in_hugo info
    cbio_not_in_hugo = hgnc.merge(cbio_genes, on=['Hugo_Symbol'], how='right')
    cbio_not_in_hugo = cbio_not_in_hugo[cbio_not_in_hugo['gene_id'].isna()]
    cbio_not_in_hugo = cbio_not_in_hugo[['Hugo_Symbol', 'Entrez_Gene_Id']]
    cbio_not_in_hugo.to_csv('counts/cbio_not_in_hugo.txt', index=None, sep='\t')

    # Generate hugo_not_in_gtf info
    hugo_not_in_gtf = final_genes[['sample', 'gene_id']].merge(hgnc, on=['gene_id'], how='right')
    hugo_not_in_gtf = hugo_not_in_gtf[hugo_not_in_gtf['sample'].isna()]
    hugo_not_in_gtf.drop(['sample'], axis=1, inplace=True)
    hugo_not_in_gtf.to_csv('counts/hugo_not_in_gtf.txt', index=None, sep='\t')
    gtf_counts = hugo_not_in_gtf.groupby(['Locus group', 'Locus type']).size()
    gtf_counts = gtf_counts.reset_index()
    gtf_counts.rename(columns={0: 'Count'}, inplace=True)
    gtf_counts.to_csv('counts/hugo_not_in_gtf_counts.txt', index=None, sep='\t')
    
    final_genes = final_genes.merge(cbio_genes, on=['Hugo_Symbol'], how='left')
    
    # Manual additions as per request
    final_genes.loc[final_genes['gene_id'] == 'ENSG00000133706', 'Hugo_Symbol'] = 'LARS'
    final_genes.loc[final_genes['gene_id'] == 'ENSG00000237452', 'Hugo_Symbol'] = 'BHMG1'

    return final_genes


def process_sample(sample, sample_data):
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


def generate_gistic_outputs(gistic_data, hdel_data, temp_dir, hgnc_file):
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
    gistic_matrix.to_csv(temp_dir + 'remixt_CNA.txt', index=None, sep='\t')


def generate_seg_outputs(aggregated_cn_data, temp_dir, stats_data):
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
        aggregated_cn_data[sample].to_csv(temp_dir + sample + '_log_seg_data.seg', index=None, sep='\t')


def clean_up_stats(stats_data):
    stats_data = pd.DataFrame(stats_data)
    stats_data['tumour_proportion'] = 1. - stats_data['normal_proportion']

    stats_data = stats_data[[
        'sample', 'ploidy', 'proportion_divergent',
        'tumour_proportion', 'proportion_divergent', 'elbo']].sort_values('sample')

    return stats_data