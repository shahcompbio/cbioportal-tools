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
import logging
import numpy as np
import os
import pandas as pd
import random
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
    ann_mapping = {
    # 'Missense_Mutation' if LOW, 'Unknown' if MODIFIER
    # Reference:
    # http://snpeff.sourceforge.net/SnpEff_manual.html ->
    # https://github.com/cbare/vcf2maf/blob/master/vcf2maf/vcf2maf.py
    # 'coding_sequence_variant': 'Missense_Mutation',
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
    # TODO: INS or DEL. Throw warning here if REF != ALT
    # (it probably shouldn't outside of strelka INDEL SNVs)
    # if REF > ALT it's a DEL, if REF < ALT it's an INS
    # DEL = Frame_Shift_Del here, INS = Frame_Shift_Ins here
    # 'frameshift_variant': 'Unknown',
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
    '3_prime_UTR_truncation': '3\'UTR', # &exon_loss
    '5_prime_UTR_variant': '5\'UTR',  
    '5_prime_UTR_truncation': '5 \'UTR',    # &exon_loss_variant
    'sequence_feature': 'Unknown'   # &exon_loss_variant
    }

    output_filename = os.path.splitext(input_file)[0].split('/')[-1] + '.maf'
    output_file = open(output_dir + output_filename, 'w+')

    output_header='#version 0.1 (cBioPortal minimal MAF format)\nHugo_Symbol\tEntrez_Gene_Id\tTumor_Sample_Barcode\tVariant_Classification\tHGVSp_Short\n'
    output_file.write(output_header)

    hugo_entrez_mapping = {}
    df = pd.read_csv(hgnc_file, delimiter='\t', dtype={'NCBI Gene ID': str})
    df = df.replace(np.nan, '')
    for hugo, entrez in zip(df['Approved symbol'], df['NCBI Gene ID']):
        hugo_entrez_mapping[hugo] = entrez

    vcf_reader = vcf.Reader(open(input_file, 'r'))
    for record in vcf_reader:
        ref = record.REF
        alt = record.ALT  # list of length 1, take length of first item
        
        record_annotations = [item.split('|')[1].split('&')[0] for item in record.INFO['ANN']]
        if 'missense_variant' in record_annotations:
            ann_index = random.choice([i for i, e in enumerate(record_annotations) if e == 'missense_variant'])
            variant_classification = 'Missense_Mutation'
            
        else:
            record_impacts = [item.split('|')[2] for item in record.INFO['ANN']]
            
            if 'HIGH' in record_impacts:
                impact = 'HIGH'
            elif 'MODERATE' in record_impacts:
                impact = 'MODERATE'
            elif 'LOW' in record_impacts:
                impact = 'LOW'
            elif 'MODIFIER' in record_impacts:
                impact = 'MODIFIER'
            else:
                logging.warning(f'{record} has no annotation impacts. Skipping record.')
                continue

            ann_index = random.choice([i for i, e in enumerate(record_impacts) if e == impact])
            annotation = record.INFO['ANN'][ann_index].split('|')[1].split('&')[0]

            if annotation == 'coding_sequence_variant':
                ann_impact = record.INFO['ANN'][ann_index].split('|')[2]
                
                if ann_impact == 'LOW':
                    variant_classification = 'Missense_Mutation'
                elif ann_impact == 'MODIFIER':
                    variant_classification = 'Unknown'
                else:
                    logging.warning(f'A variant_classification for {annotation} with impact {ann_impact} was not found.')
                    variant_classification = 'Unknown'

            elif annotation == 'frameshift_variant':
                if len(ref) != len(alt[0]):
                    logging.warning(f'{record} has REF and ALT with unequal lengths')
                    if len(ref) > len(alt[0]):
                        variant_classification = 'Frame_Shift_Del'
                    elif len(ref) < len(alt[0]):
                        variant_classification = 'Frame_Shift_Ins'
                else:
                    variant_classification = 'Unknown'
            else:
                if annotation in ann_mapping:
                    variant_classification = ann_mapping[annotation]
                else:
                    logging.warning(f'A variant_classification for {annotation} was not found.')
                    variant_classification = 'Unknown'

        record_ann = record.INFO['ANN'][ann_index].split('|')
        hugo_symbol = record_ann[3]
        
        if hugo_symbol in hugo_entrez_mapping:
            entrez_gene_id = hugo_entrez_mapping[hugo_symbol]
        else:
            entrez_gene_id = ''

        hgvsp_short = record_ann[10]
            
        output_file.write(hugo_symbol + '\t' + entrez_gene_id + '\t' + sample_id + '\t' + variant_classification +'\t' + hgvsp_short + '\n')


@click.command()
@click.argument('input_file')
@click.argument('sample_id')
@click.argument('hgnc_file')
@click.option('--output_dir', default='')
def main(input_file, sample_id, hgnc_file, output_dir):
    convert(input_file, sample_id, hgnc_file, output_dir)


if __name__ == '__main__':
    main()
