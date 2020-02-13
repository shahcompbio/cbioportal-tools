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


@click.command()
@click.argument('input_yaml')
@click.argument('path_to_output_study')
@click.option('--temp_dir', default='')
@click.option('--path_to_external_study')
def main(input_yaml, path_to_output_study, temp_dir, path_to_external_study):
    with open(input_yaml) as file:
        yaml_file = yaml.full_load(file)
        
        path_to_mapping_file = yaml_file['id_mapping']
        path_to_gtf_file = yaml_file['gtf']
        
        for patient, doc in yaml_file['patients'].items():
            for sample, doc in doc.items():
                filter(sample, doc['museq_vcf'], doc['strelka_vcf'], temp_dir)


if __name__ == '__main__':
    main()
