import click

from generate_outputs import extract, transform, load


@click.command()
@click.argument('gtf')
@click.argument('hgnc')
@click.argument('igv_segs')
@click.argument('titan_segs')
@click.argument('sample_id')
@click.option('--output_dir', default='test/test_output/')
def main(gtf, hgnc, igv_segs, titan_segs, sample_id, output_dir):
    extract(gtf, hgnc, igv_segs, titan_segs)
    gene_dict, seg_dict = transform()
    load(gene_dict, seg_dict, sample_id, output_dir)


if __name__ == '__main__':
    main()
