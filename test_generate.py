import click
import filecmp
import shutil

from generate_outputs import extract, transform, load


# check outputs and intermediate file against baselines
@click.command()
@click.argument('gtf')
@click.argument('hgnc')
@click.argument('igv_segs')
@click.argument('titan_segs')
@click.argument('sample_id')
@click.option('--output_dir', default='test/generate_outputs/test_output/')
@click.option('--show_missing_hugo/--no_missing_hugo', default=False)
@click.option('--show_missing_entrez/--no_missing_entrez', default=False)
@click.option('--show_missing_both/--no_missing_both', default=False)
def main(gtf, hgnc, igv_segs, titan_segs, sample_id, output_dir, show_missing_hugo, show_missing_entrez, show_missing_both):
    extracted_file = extract(gtf, hgnc, igv_segs, titan_segs)
    gene_dict, seg_dict = transform(extracted_file, show_missing_hugo, show_missing_entrez, show_missing_both)
    load(gene_dict, seg_dict, sample_id, output_dir)

    if filecmp.cmp('test/generate_outputs/output_baseline/gistic_gene_data.txt', output_dir + 'gistic_gene_data.txt'):
        print('Gistic gene data output matches baseline.')
    else:
        print('Gistic gene data output does not match baseline!')
    
    if filecmp.cmp('test/generate_outputs/output_baseline/log_seg_data.txt', output_dir + 'log_seg_data.txt'):
        print('Log segment data output matches baseline.')
    else:
        print('Log segment data output does not match baseline!')
    
    if filecmp.cmp('test/generate_outputs/output_baseline/integer_gene_data.txt', output_dir + 'integer_gene_data.txt'):
        print('Integer gene data output matches baseline.')
    else:
        print('Integer gene data output does not match baseline!')
    
    if filecmp.cmp('test/generate_outputs/output_baseline/integer_seg_data.txt', output_dir + 'integer_seg_data.txt'):
        print('Integer segment data output matches baseline.')
    else:
        print('Integer segment data output does not match baseline!')


if __name__ == '__main__':
    main()
