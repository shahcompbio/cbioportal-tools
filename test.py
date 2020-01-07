import click
import filecmp

from generate_outputs import extract, transform, load


@click.command()
@click.argument('gtf')
@click.argument('hgnc')
@click.argument('igv_segs')
@click.argument('titan_segs')
@click.argument('sample_id')
@click.option('--output_dir', default='test/test_output/')
def main(gtf, hgnc, igv_segs, titan_segs, sample_id, output_dir):
    extracted_file = extract(gtf, hgnc, igv_segs, titan_segs)
    gene_dict, seg_dict = transform(extracted_file)
    load(gene_dict, seg_dict, sample_id, output_dir)

    print('Outcomes of comparing test results to baselines:')

    if filecmp.cmp('test/output_baseline/gistic_gene_data.txt', output_dir + 'gistic_gene_data.txt'):
        print('Gistic gene data output matches baseline.')
    else:
        print('Gistic gene data output does not match baseline!')
    
    if filecmp.cmp('test/output_baseline/gistic_seg_data.txt', output_dir + 'gistic_seg_data.txt'):
        print('Gistic segment data output matches baseline.')
    else:
        print('Gistic segment data output does not match baseline!')
    
    if filecmp.cmp('test/output_baseline/integer_gene_data.txt', output_dir + 'integer_gene_data.txt'):
        print('Integer gene data output matches baseline.')
    else:
        print('Integer gene data output does not match baseline!')
    
    if filecmp.cmp('test/output_baseline/integer_seg_data.txt', output_dir + 'integer_seg_data.txt'):
        print('Integer segment data output matches baseline.')
    else:
        print('Integer segment data output does not match baseline!')


if __name__ == '__main__':
    main()
