import click
import glob
import pandas as pd

# merge multiple gistic OR integer gene data text files
def merge_gene_outputs(input_dir, output_dir):
    files_to_merge = glob.glob(input_dir + '*.txt')
    files_with_issues = []
    dfs_to_merge = []
    
    for file in files_to_merge:
        data_frame = pd.read_csv(file, delimiter='\t', dtype={'entrez_id': str}).astype(str)
        
        if list(data_frame)[0] != 'entrez_id' or list(data_frame)[1] != 'hugo_symbol':
            files_with_issues.append(file)
            continue
        
        dfs_to_merge.append(data_frame)

    if files_with_issues:
        print('The following list contains files with issues:')
        print(files_with_issues)
        print('Perhaps they are segment data, or have malformed headers?')
        print('Please fix or remove them and re-run the script.')
        return

    merged_file = dfs_to_merge.pop()
    
    while dfs_to_merge:
        merged_file = pd.merge(merged_file, dfs_to_merge.pop(), on=['hugo_symbol', 'entrez_id'], how='outer')

    merged_file.to_csv(output_dir + 'merged.txt', index=None, sep='\t')


@click.command()
@click.option('--input_dir', default='')
@click.option('--output_dir', default='')
def main(input_dir, output_dir):
    merge_gene_outputs(input_dir, output_dir)


if __name__ == '__main__':
    main()
