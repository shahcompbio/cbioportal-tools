import click
import filecmp

from merge_gene_outputs import merge_gene_outputs


@click.command()
@click.option('--input_dir', default='')
@click.option('--output_dir', default='')
def main(input_dir, output_dir):
    merge_gene_outputs(input_dir, output_dir)


if __name__ == '__main__':
    main()