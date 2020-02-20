'''
    for each section in yaml, perform:
        filtering ✔
        run cn data output generator ✔
        conversion to maf

    afterwards:
        merge mafs ✔
        merge gistic gene results ✔
        merge log seg data results ✔
        create cBio stidy meta files/connect to data files properly

    merge wgs with impact data?
    idea:   take two existing studies and compare files to see if any
            any overlap, take new description, name etc. as user
            input and merge exisitng non-overlapping data.
    '''

import click
import gzip
import yaml

from convert_vcf_to_maf import convert as convert_vcf_to_maf
from generate_outputs import extract, transform, load
from merge_outputs import merge_all_data as merge_outputs


def merge_studies(path_to_output_study, path_to_external_study, output_dir)
    pass


def create_study(work_dir, path_to_output_study):
    pass


def generate_outputs(gtf_file, hgnc_file, titan_igv, titan_segs, sample_id, output_dir): 
    extracted_file = extract(gtf_file, hgnc_file, titan_igv, titan_segs)
    gene_dict, seg_dict = transform(extracted_file, show_missing_hugo=False, show_missing_entrez=False, show_missing_both=False)
    load(gene_dict, seg_dict, sample_id, output_dir, output_gistic_gene=True, output_integer_gene=False, output_log_seg=True, output_integer_seg=False)


def filter(sample_id, museq_vcf, strelka_vcf, work_dir):
    '''
    original code by Diljot Grewal

    required mofifications (done): take museq and strelka directly as input, output to temp_dir
    strip sample id from inputs and append to output filtered file
    '''

    museq_filtered = work_dir + '{}_museq_filtered.vcf.gz'.format(sample_id)

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


@click.command()
@click.argument('input_yaml')
@click.argument('path_to_output_study')
@click.option('--temp_dir', default='')
@click.option('--path_to_external_study')
@click.option('--path_to_merged_study')
def main(input_yaml, path_to_output_study, temp_dir, path_to_external_study):
    with open(input_yaml) as file:
        yaml_file = yaml.full_load(file)
        hgnc_file = yaml_file['id_mapping']
        gtf_file = yaml_file['gtf']
        
        for patient, doc in yaml_file['patients'].items():
            for sample, doc in doc.items():
                museq_filtered = filter(sample, doc['museq_vcf'], doc['strelka_vcf'], temp_dir)
                
                with gzip.open(museq_filtered, 'rt') as museq_vcf:
                    # TODO: write converter
                    convert_vcf_to_maf(museq_vcf, sample, hgnc_file, temp_dir)
                
                with gzip.open(doc['titan_igv'], 'rt') as titan_igv, gzip.open(doc['titan_segs'], 'rt') as titan_segs:
                    generate_outputs(gtf_file, hgnc_file, titan_igv, titan_segs, sample, temp_dir)
                
    merge_outputs(temp_dir, temp_dir)
    create_study(temp_dir, path_to_output_study)

    if path_to_external_study:
        merge_studies(path_to_output_study, path_to_external_study, path_to_merged_study)


if __name__ == '__main__':
    main()
