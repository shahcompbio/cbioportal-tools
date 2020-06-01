import click
import os
import subprocess


def convert(input_file, sample_id, dataset_id, output_dir):
    assert input_file.endswith('.gz')

    vcf_filename = os.path.join(output_dir, f'{dataset_id}.vcf')
    maf_filename = os.path.join(output_dir, f'{dataset_id}-generated.maf')

    with open(vcf_filename, 'w') as f:
        subprocess.check_call(['gunzip', '-c', input_file], stdout=f)

    cmd = (f'/opt/local/singularity/3.3.0/bin/singularity run --bind /juno/work/shah/svatrt/vcf2maf:/vcf2maf \
            docker://wgspipeline/vcf2maf:v0.0.1 \
            vcf2maf.pl \
            --input-vcf {vcf_filename} \
            --output-maf {maf_filename} \
            --vep-path /usr/local/bin \
            --ref-fasta /vcf2maf/cache/homo_sapiens/99_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
            --filter-vcf /vcf2maf/cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
            --vep-data /vcf2maf/cache/ \
            --tumor-id {sample_id}')

    subprocess.check_call(cmd.split())


@click.command()
@click.argument('input_file')
@click.argument('sample_id')
@click.argument('dataset_id')
@click.option('--output_dir', default='')
def main(input_file, sample_id, dataset_id, output_dir):
    convert(input_file, sample_id, dataset_id, output_dir)


if __name__ == '__main__':
    main()
