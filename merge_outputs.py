import click
import glob
import numpy as np
import os
import pandas as pd

from collections import Counter

# merge multiple gistic OR integer gene data text files
def merge_gistic_gene_data(input_dir, output_dir):
    files_to_merge = [fn for fn in glob.glob(input_dir + '*.txt') if not os.path.basename(fn).startswith('merged')]
    files_to_merge = sorted(files_to_merge)
    files_with_issues = []
    dfs_to_merge = []

    for file in files_to_merge:
        data_frame = pd.read_csv(file, delimiter='\t', dtype=str)
        
        if list(data_frame)[0] != 'Hugo_Symbol' or list(data_frame)[1] != 'Entrez_Gene_Id':
            files_with_issues.append(file)
            continue
        
        dfs_to_merge.append(data_frame)

    if files_with_issues:
        print('The following list contains gistic gene data files that could not be merged:')
        print(files_with_issues)
        print('Please fix or remove them and re-run the merge script.')
        return

    merged_file = dfs_to_merge.pop(0)
    
    while dfs_to_merge:
        merged_file = pd.merge(merged_file, dfs_to_merge.pop(), on=['Hugo_Symbol', 'Entrez_Gene_Id'], how='outer')

    merged_file = merged_file.replace(np.nan, 'NA')
    merged_file = merged_file.replace('nan', '')
    merged_file.to_csv(output_dir + 'merged.txt', index=None, sep='\t')


def merge_log_seg_data(input_dir, output_dir):
    # skip 1 line in each file after the first
    files_to_merge = [fn for fn in glob.glob(input_dir + '*.seg') if not os.path.basename(fn).startswith('merged')]
    files_to_merge = sorted(files_to_merge)
    
    with open(output_dir + 'merged.seg', 'w+') as outfile:
        with open(files_to_merge.pop(0)) as infile:
            outfile.write(infile.read())
        for file in files_to_merge:
            with open(file) as infile:
                next(infile)
                outfile.write(infile.read())


def merge_maf_data(input_dir, output_dir):
    # skip 2 lines in each file after the first
    files_to_merge = [fn for fn in glob.glob(input_dir + '*.maf') if not os.path.basename(fn).startswith('merged')]
    files_to_merge = sorted(files_to_merge)
    
    with open(output_dir + 'merged.maf', 'w+') as outfile:
        with open(files_to_merge.pop(0)) as infile:
            outfile.write(infile.read())
        for file in files_to_merge:
            with open(file) as infile:
                for line in infile.readlines()[2:]:
                    outfile.write(line)


def merge_all_data(input_dir, output_dir):
    merge_gistic_gene_data(input_dir, output_dir)
    merge_log_seg_data(input_dir, output_dir)
    merge_maf_data(input_dir, output_dir)


@click.command()
@click.option('--file_types', '-ft', type=click.Choice(['all', 'gistic_gene', 'log_seg', 'maf']), multiple=True, default=['all'])
@click.option('--input_dir', default='')
@click.option('--output_dir', default='')
def main(file_types, input_dir, output_dir):
    if 'all' in file_types or Counter(file_types) == Counter(['gistic_gene', 'log_seg', 'maf']):
        merge_all_data(input_dir, output_dir)
        return
    
    if 'gistic_gene' in file_types:
        merge_gistic_gene_data(input_dir, output_dir)

    if 'log_seg' in file_types:
        merge_log_seg_data(input_dir, output_dir)

    if 'maf' in file_types:
        merge_maf_data(input_dir, output_dir)


if __name__ == '__main__':
    main()
