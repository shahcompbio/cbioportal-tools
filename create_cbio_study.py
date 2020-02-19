'''
    for each section in yaml, perform:
        filtering
        run cn data output generator
        conversion to maf

    afterwards:
        merge mafs
        merge gistic gene results
        merge log seg data results
        create cBio stidy meta files/connect to data files properly

    merge wgs with impact data?
    '''

import click
import gzip
import yaml

from convert_vcf_to_maf import convert
from generate_outputs import extract, transform, load
from merge_outputs import merge_gistic_gene_data, merge_log_seg_data, merge_maf_data, merge_all_data


def filter(sample_id, museq_vcf, strelka_vcf, temp_dir):
    '''
    original code by Diljot Grewal

    required mofifications (done): take museq and strelka directly as input, output to temp_dir
    strip sample id from inputs and append to output filtered file
    '''

    museq_filtered = temp_dir + '{}_museq_filtered.vcf.gz'.format(sample_id)

    strelka_ref = set()

    with gzip.open(strelka_vcf, 'rt') as strelka_data:
        for line in strelka_data:
            if line.startswith('#'):
                continue
            
            line = line.strip().split()
            
            chrom = line[0]
            pos = line[1]
            
            strelka_ref.add((chrom, pos))

    with gzip.open(museq_vcf, 'rt') as museq_data, gzip.open(museq_filtered, 'wt') as museqout:
        for line in museq_data:
            if line.startswith('#'):
                museqout.write(line)
                continue
            
            line = line.strip().split()
            chrom = line[0]
            pos = line[1]
            
            if ((chrom, pos))  not in strelka_ref:
                continue
            
            pr = line[7].split(';')[0].split('=')[1]
            if float(pr) < 0.85:
                continue
            
            outstr = '\t'.join(line)+'\n'
            museqout.write(outstr)

    return museq_filtered


cli = click.Group()

@cli.command()
@click.argument('input_file')
@click.argument('sample_id')
@click.argument('hgnc_file')
@click.option('--output_dir', default='')
def convert_vcf_to_maf(input_file, sample_id, hgnc_file, output_dir):
    convert(input_file, sample_id, hgnc_file, output_dir)


@cli.command()
@click.argument('gtf_file')
@click.argument('hgnc_file')
@click.argument('titan_igv')
@click.argument('titan_segs')
@click.argument('sample_id')
@click.option('--output_dir', default='')
def main(gtf_file, hgnc_file, titan_igv, titan_segs, sample_id, output_dir): 
    extracted_file = extract(gtf_file, hgnc_file, titan_igv, titan_segs)
    gene_dict, seg_dict = transform(extracted_file, show_missing_hugo=False, show_missing_entrez=False, show_missing_both=False)
    load(gene_dict, seg_dict, sample_id, output_dir, output_gistic_gene=True, output_integer_gene=False, output_log_seg=True, output_integer_seg=False)


@cli.command()
@click.option('--input_dir', default='')
@click.option('--output_dir', default='')
def merge_outputs(file_types, input_dir, output_dir):
    if 'all' in file_types or Counter(file_types) == Counter(['gistic_gene', 'log_seg', 'maf']):
        merge_all_data(input_dir, output_dir)
        return
    
    if 'gistic_gene' in file_types:
        merge_gistic_gene_data(input_dir, output_dir)

    if 'log_seg' in file_types:
        merge_log_seg_data(input_dir, output_dir)

    if 'maf' in file_types:
        merge_maf_data(input_dir, output_dir)


@cli.command()
@click.argument('input_yaml')
@click.argument('path_to_output_study')
@click.option('--temp_dir', default='')
@click.option('--path_to_external_study')
@click.pass_context
def main(ctx, input_yaml, path_to_output_study, temp_dir, path_to_external_study):
    with open(input_yaml) as file:
        yaml_file = yaml.full_load(file)
        hgnc_file = yaml_file['id_mapping']
        gtf_file = yaml_file['gtf']
        
        for patient, doc in yaml_file['patients'].items():
            for sample, doc in doc.items():
                museq_filtered = filter(sample, doc['museq_vcf'], doc['strelka_vcf'], temp_dir)
                
                with gzip.open(museq_filtered, 'rt') as museq_vcf:
                    # TODO: write converter
                    ctx.invoke(convert_vcf_to_maf, input_file=museq_filtered, sample_id=sample, hgnc_file=hgnc_file, output_dir=temp_dir)
                
                with gzip.open(doc['titan_igv'], 'rt') as titan_igv, gzip.open(doc['titan_segs'], 'rt') as titan_segs:
                    ctx.invoke(generate_outputs, gtf_file=gtf_file, hgnc_file=hgnc_file, titan_igv=titan_igv, titan_segs=titan_segs, sample_id=sample, output_dir=temp_dir)
                
                ctx.invoke(merge_outputs, input_dir=temp_dir, output_dir=temp_dir)


if __name__ == '__main__':
    main()
