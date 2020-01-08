import click
import filecmp

from merge_gene_outputs import merge_gene_outputs


@click.command()
@click.option('--input_dir', default='test/merge_gene_outputs/test_input/')
@click.option('--output_dir', default='test/merge_gene_outputs/test_output/')
def main(input_dir, output_dir):
    merge_gene_outputs(input_dir, output_dir)

    if filecmp.cmp('test/merge_gene_outputs/output_baseline/merged.txt', output_dir + 'merged.txt'):
        print('Merged file matches baseline.')
    else:
        print('Merged file does not matche baseline!')


if __name__ == '__main__':
    main()
