import pandas as pd
import numpy as np


autosomes = [str(a) for a in range(1, 23)]


def aggregate_adjacent(cnv, value_cols=(), stable_cols=(), length_normalized_cols=(), summed_cols=()):
    """ Aggregate adjacent segments with similar copy number state.

    see: https://github.com/amcpherson/remixt/blob/master/remixt/segalg.py

    Args:
        cnv (pandas.DataFrame): copy number table

    KwArgs:
        value_cols (list): list of columns to compare for equivalent copy number state
        stable_cols (list): columns for which values are the same between equivalent states
        length_normalized_cols (list): columns that are width normalized for equivalent states

    Returns:
        pandas.DataFrame: copy number with adjacent segments aggregated
    """

    # Group segments with same state
    cnv = cnv.sort_values(['chr'] + value_cols)
    cnv['chromosome_index'] = np.searchsorted(np.unique(cnv['chr']), cnv['chr'])
    cnv['diff'] = cnv[['chromosome_index'] + value_cols].diff().abs().sum(axis=1)
    cnv['is_diff'] = (cnv['diff'] != 0)
    cnv['cn_group'] = cnv['is_diff'].cumsum()

    def agg_segments(df):
        a = df[stable_cols].iloc[0]

        a['chr'] = df['chr'].min()
        a['start'] = df['start'].min()
        a['end'] = df['end'].max()
        a['width'] = df['width'].sum()

        for col in length_normalized_cols:
            a[col] = (df[col] * df['width']).sum() / (df['width'].sum() + 1e-16)

        for col in summed_cols:
            a[col] = df[col].sum()

        return a

    aggregated = cnv.groupby('cn_group').apply(agg_segments)

    for col in aggregated:
        aggregated[col] = aggregated[col].astype(cnv[col].dtype)

    return aggregated


def calculate_gene_copy(cnv, genes):
    """ Calculate the copy number segments overlapping each gene

    Args:
        cnv (pandas.DataFrame): copy number table
        genes (pandas.DataFrame): gene table

    Returns:
        pandas.DataFrame: segment copy number for each gene

    The input copy number table is assumed to have columns: 
        'chr', 'start', 'end', 'width', 'copy', 'reads', 'state'

    The input genes table is assumed to have columns:
        'chr', 'gene_start', 'gene_end', 'gene_id'

    The output segment copy number table should have columns:
        'gene_id', 'chr', 'gene_start', 'gene_end', 'start', 'end',
        'width', 'copy', 'reads', 'state'
    where each entry in the output table represents an overlap between
    a gene and a segment.

    """

    data = []

    for chr in cnv['chr'].unique():
        chr_cnv = cnv[cnv['chr'] == chr]
        chr_genes = genes[genes['chr'] == chr]

        # Iterate through segments, calculate overlapping genes
        for idx, row in chr_cnv.iterrows():
            
            # Subset overlapping genes
            overlapping_genes = chr_genes[~((chr_genes['gene_end'] < row['start']) | (chr_genes['gene_start'] > row['end']))]

            # Add cnv columns
            overlapping_genes

            data.append(overlapping_genes)

    data = pd.concat(data, ignore_index=True)

    return data


def test_calculate_gene_copy():
    # Create some cnv data
    cnv = pd.DataFrame({})

    # Create some gene data
    gene = pd.DataFrame({})

    # Create known overlap
    overlapping = pd.DataFrame({})

    assert overlapping == calculate_gene_copy(cnv, genes)


def read_copy_data(bins_filename, filter_normal=False):
    """ Read hmmcopy data, filter normal cells and aggregate into segments
    """
    data = pd.read_csv(bins_filename)

    # Filter normal cells that are approximately diploid
    if filter_normal:
        cell_stats = (
            data[data['chr'].isin(autosomes)]
            .groupby('cell_id')['state']
            .agg(['mean', 'std'])
            .reset_index())

        normal_cells = (
            cell_stats
            .query('1.95 < mean < 2.05')
            .query('std < 0.01'))

        data = data.merge(normal_cells[['cell_id']])

    # Aggregate cell copy number
    data = (
        data
        .groupby(['chr', 'start', 'end', 'width'])
        .agg({'state': 'median', 'copy': np.nanmean, 'reads': 'sum'})
        .reset_index())

    # Aggregate cell copy number
    data = aggregate_adjacent(
        data,
        value_cols=['state'],
        stable_cols=['state'],
        length_normalized_cols=['copy'],
        summed_cols=['reads'],
    )

    return data
