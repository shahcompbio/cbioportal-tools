'''
References:
1. http://snpeff.sourceforge.net/SnpEff_manual.html
2. https://github.com/cbare/vcf2maf/blob/master/vcf2maf/vcf2maf.py
3. https://svn.bcgsc.ca/bitbucket/projects/KRONOS/repos/convert_vcf_to_maf/browse/component_seed/vcf2maf/vcf2maf-master/vcf2maf_museq.pl

If two or more annotations in .vcf Annotation column, take first one
(maybe related to "most deleterious" as shown here:
http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf)

Annotation_Impact (impact prediction):
HIGH > MODERATE > LOW > MODIFIER

If Annotation column in .vcf record INFO column is 'missense_variant'
treat as highest rank (above HIGH)
(https://github.com/shahcompbio/scgenome/pull/16)

If Annotation column is 'coding_sequence_variant':
if Annotation_Impact is 'LOW', map to 'Missense_Mutation',
if Annotation_Impact is 'MODIFIER', map to 'Unknown'
Reference: link #2

If Annotation column is 'frameshift_variant':
inspect REF and ALT .vcf columns and ompare number of characters
if REF > ALT map to Frame_Shift_Del,
if REF < ALT map to Frame_Shift_Ins,
otherwise map to 'Unknown'
REF and ALT length should be the same outside of strelka INDEL SNVs
Reference: Diljot Grewal

If multiple annotations tie for highest rank, pick one at random
'''

import click
import csv
import glob
import logging
import numpy as np
import pandas as pd
import random
import vcf


def convert(input_file, sample_id, hgnc_file, output_dir):
    '''
    required columns:
        'Hugo_Symbol' (Gene_Name)
        'Entrez_Gene_Id' (map from custom.txt)
        'Tumor_Sample_Barcode' (user input)
        'Variant_Classification' (from Annotation and ann_mapping dict)
        'HGVSp_Short' (HGVS.p)
    '''
    
    ann_mapping = {
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
    'structural_interaction_variant': 'Unknown', # my own conclusion
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
    '5_prime_UTR_truncation': '5 \'UTR', # &exon_loss_variant
    'sequence_feature': 'Unknown' # &exon_loss_variant
    }

    output_filename = '{}_strelka_indel_annotated.maf'.format(sample_id)
    output_file = open(output_dir + output_filename, 'w+')

    output_header='#version 0.1 (cBioPortal minimal MAF format)\nHugo_Symbol\tEntrez_Gene_Id\tVariant_Classification\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode\tHGVSp_Short\n'
    output_file.write(output_header)

    hugo_entrez_mapping = {}
    df = pd.read_csv(hgnc_file, delimiter='\t', dtype={'NCBI Gene ID': str})
    df = df.replace(np.nan, '')
    for hugo, entrez in zip(df['Approved symbol'], df['NCBI Gene ID']):
        hugo_entrez_mapping[hugo] = entrez
    
    print(f'Running conversion for sample {sample_id}.')

    vcf_reader = vcf.Reader(filename=input_file)
    for record in vcf_reader:
        ref = record.REF
        alt = ''.join(map(str, record.ALT))
        
        record_annotations = [ann.split('|')[1].split('&')[0] for ann in record.INFO['ANN']]
        if 'missense_variant' in record_annotations:
            ann_index = random.choice([i for i, e in enumerate(record_annotations) if e == 'missense_variant'])
            variant_classification = 'Missense_Mutation'
            
        else:
            record_impacts = [ann.split('|')[2] for ann in record.INFO['ANN']]
            
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
                if len(ref) != len(alt):
                    logging.warning(f'{record} has REF and ALT with unequal lengths!')
                    if len(ref) > len(alt):
                        variant_classification = 'Frame_Shift_Del'
                    elif len(ref) < len(alt):
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
        if record_ann[3] == '':
            hugo_symbol = 'Unknown'
        else:
            hugo_symbol = record_ann[3]
        
        if hugo_symbol in hugo_entrez_mapping:
            if hugo_entrez_mapping[hugo_symbol] == '':
                entrez_gene_id = '0'
            else:
                entrez_gene_id = hugo_entrez_mapping[hugo_symbol]
        else:
            entrez_gene_id = '0'

        hgvsp_short = record_ann[10]
            
        output_file.write(hugo_symbol + '\t' + entrez_gene_id + '\t' + variant_classification + '\t' + ref + '\t' + ref + '\t' + alt + '\t' + sample_id +'\t' + hgvsp_short + '\n')
    
    print(f'Conversion finished for sample {sample_id}.')


@click.command()
@click.argument('input_file')
@click.argument('sample_id')
@click.argument('hgnc_file')
@click.option('--output_dir', default='')
def main(input_file, sample_id, hgnc_file, output_dir):
    convert(input_file, sample_id, hgnc_file, output_dir)


if __name__ == '__main__':
    main()
