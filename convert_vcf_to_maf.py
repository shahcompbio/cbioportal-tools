'''
http://snpeff.sourceforge.net/SnpEff_manual.html
https://github.com/cbare/vcf2maf/blob/master/vcf2maf/vcf2maf.py
https://svn.bcgsc.ca/bitbucket/projects/KRONOS/repos/convert_vcf_to_maf/browse/component_seed/vcf2maf/vcf2maf-master/vcf2maf_museq.pl
if two or more annotations in .vcf Annotation column, take first one (maybe related to "most deleterious" as shown here: http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf) 
'''

import click
import csv
import glob
import pandas as pd

@click.command()
@click.option('--input_dir', default='')
@click.option('--output_dir', default='')
def main(input_dir, output_dir):
    required_cols = [
    'Hugo_Symbol', # Gene_Name
    'Entrez_Gene_Id', # map from custom.txt
    'Tumor_Sample_Barcode', # user input
    'Variant_Classification',  # need to develop mapping (Annotation)
    'HGVSp_Short', # HGVS.p
    ]
    
    vcf_to_maf = {
    'coding_sequence_variant': 'Missense_Mutation', # LOW, 'Unknown' if MODIFIER -> http://snpeff.sourceforge.net/SnpEff_manual.html
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
    'frameshift_variant': 'Unknown',    # TODO: INS or DEL. In strelka INDEL SNVs, throw warning or exception here if this occurs
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
    '3_prime_UTR_truncation&exon_loss': '3\'UTR',   # &exon_loss
    '5_prime_UTR_variant': '5\'UTR',  
    '5_prime_UTR_truncation&exon_loss_variant': '5 \'UTR', # &exon_loss_variant
    'sequence_feature&exon_loss_variant': 'Unknown'    # &exon_loss_variant
    }

    # TODO: if not in the mapping, throw some kind of exception/warning

    maf_choices = [
    'Frame_Shift_Del',
    'Frame_Shift_Ins',
    'In_Frame_Del',
    'In_Frame_Ins',
    'Missense_Mutation',
    'Nonsense_Mutation',
    'Silent',
    'Splice_Site',
    'Translation_Start_Site',
    'Nonstop_Mutation',
    '3\'UTR',
    '3\'Flank',
    '5\'UTR',
    '5\'Flank',
    'IGR',
    'Intron',
    'RNA',
    'Targeted_Region',
    'De_novo_Start_InFrame',
    'De_novo_Start_OutOfFrame',
    'Splice_Region',
    'Unknown'
    ]

    # 'MODIFER' at column 92
    open_file = open(input_dir + 'data_mutations_extended.maf', 'r')
    next(open_file)
    next(open_file)
    read = csv.reader(open_file, delimiter='\t')

    df = pd.read_csv(input_dir + 'data_mutations_extended.maf', delimiter='\t', skiprows=1)
    df = df[required_cols]
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', 2000)

    files_to_merge = glob.glob(input_dir + '*.vcf')
    # for file in files_to_merge:          


if __name__ == '__main__':
    main()
