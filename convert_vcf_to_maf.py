'''
http://snpeff.sourceforge.net/SnpEff_manual.html
https://github.com/cbare/vcf2maf/blob/master/vcf2maf/vcf2maf.py
https://svn.bcgsc.ca/bitbucket/projects/KRONOS/repos/convert_vcf_to_maf/browse/component_seed/vcf2maf/vcf2maf-master/vcf2maf_museq.pl

if two or more annotations in .vcf Annotation column, take first one (maybe related to "most deleterious" as shown here: http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf)

Impact prediction:
HIGH > MODERATE > LOW > MODIFIER
If Annotation column in .vcf is 'missense_variant', treat as highest rank (above HIGH) (https://github.com/shahcompbio/scgenome/pull/16)
If multiple annotations tie for highest rank, pick one at random
'''

import click
import csv
import glob
import numpy as np
import os
import pandas as pd
import vcf


def convert(input_file, sample_id, hgnc_file, output_dir):
    required_cols = [
    'Hugo_Symbol', # Gene_Name
    'Entrez_Gene_Id', # map from custom.txt
    'Tumor_Sample_Barcode', # user input
    'Variant_Classification',  # need to develop mapping (Annotation)
    'HGVSp_Short', # HGVS.p
    ]
    
    # TODO: if not in the mapping, throw warning, mark as 'Unknown'
    vcf_to_maf = {
    # LOW, 'Unknown' if MODIFIER
    # Reference: http://snpeff.sourceforge.net/SnpEff_manual.html
    'coding_sequence_variant': 'Missense_Mutation',
    'chromosome': 'Unknown',
    'duplication': 'Unknown',
    'inversion': 'Unknown',
    'inframe_insertion': 'In_Frame_Ins',
    'disruptive_inframe_insertion': 'In_Frame_Ins',
    'inframe_deletion': 'In_Frame_Del',
    'disruptive_inframe_deletion': 'In_Frame_Del',
    'downstream_gene_variant': '3\'Flank',
    'exon_variant': 'Unknown',
    'exon_loss_variant': 'In_Frame_Del',
    # TODO: INS or DEL. Throw warning here if this occurs
    # (it probably shouldn't outside of strelka INDEL SNVs)
    # if REF > ALT it's a DEL, if REF < ALT it's an INS
    # DEL = Frame_Shift_Del here, INS = Frame_Shift_Ins here
    'frameshift_variant': 'Unknown',
    'gene_variant': 'Unknown',
    'feature_ablation': 'Unknown',
    'gene_fusion': 'Unknown',
    'bidirectional_gene_fusion': 'Unknown',
    'rearranged_at_DNA_level': 'Unknown',
    'intergenic_region': 'IGR',
    'conserved_intergenic_variant': 'IGR',
    'intragenic_variant': 'IGR',
    'intron_variant': 'Intron',
    'conserved_intron_variant': 'Intron',
    'miRNA': 'RNA', # my own conclusion
    'missense_variant': 'Missense_Mutation',
    'initiator_codon_variant': 'Translation_Start_Site',
    'stop_retained_variant': 'Unknown',
    'protein_protein_contact': 'Unknown',
    'structural_interaction_variant': '',
    'rare_amino_acid_variant': 'Missense_Mutation',
    'splice_acceptor_variant': 'Splice_Site',
    'splice_donor_variant': 'Splice_Site',
    'splice_region_variant': 'Splice_Site', # my own conclusion
    'stop_lost': 'Nonstop_Mutation',
    '5_prime_UTR_premature_start_codon_gain_variant': 'Translation_Start_Site',
    'start_lost': 'Translation_Start_Site',
    'stop_gained': 'Missense_Mutation',
    'synonymous_variant': 'Silent',
    'start_retained': 'Translation_Start_Site',
    'stop_retained_variant': 'Silent',
    'transcript_variant': 'Unknown',
    'regulatory_region_variant': 'Unknown',
    'upstream_gene_variant': '5\'Flank',
    '3_prime_UTR_variant': '3\'UTR',
    '3_prime_UTR_truncation&exon_loss': '3\'UTR',
    '5_prime_UTR_variant': '5\'UTR',  
    '5_prime_UTR_truncation&exon_loss_variant': '5 \'UTR',
    'sequence_feature&exon_loss_variant': 'Unknown'
    }

    header='#version 0.1 (cBioPortal minimal MAF format)\nHugo_Symbol\tEntrez_Gene_Id\tTumor_Sample_Barcode\tVariant_Classification\tHGVSp_Short\n'
    output_filename = os.path.splitext(input_file)[0].split('/')[-1] + '.maf'

    hugo_entrez_mapping = {}
    df = pd.read_csv(hgnc_file, delimiter='\t', dtype={'NCBI Gene ID': str})
    df = df.replace(np.nan, '')
    for hugo, entrez in zip(df['Approved symbol'], df['NCBI Gene ID']):
        hugo_entrez_mapping[hugo] = entrez

    vcf_reader = vcf.Reader(open(input_file, 'r'))
    for record in vcf_reader:
        record.REF
        record.ALT  # list of length 1, take length of first item
        for item in record.INFO['ANN']:
            spl = item.split('|')
            spl[1]  # Annotation
            spl[2]  # Annotation_Impact 
            spl[3]  # Gene_Name
            spl[10] # HGVS.p, missing often in twins. SnpEFF docs?


@click.command()
@click.argument('input_file')
@click.argument('sample_id')
@click.argument('hgnc_file')
@click.option('--output_dir', default='')
def main(input_file, sample_id, hgnc_file, output_dir):
    convert(input_file, sample_id, hgnc_file, output_dir)


if __name__ == '__main__':
    main()
