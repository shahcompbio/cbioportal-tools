# original code by Diljot Grewal

import gzip
import sys

sampleid = sys.argv[1]

museq = '{}_museq_paired_annotated.vcf.gz'.format(sampleid)
strelka = '{}_strelka_snv_annotated.vcf.gz'.format(sampleid)
museq_filtered = '{}_museq_filtered.vcf.gz'.format(sampleid)

strelka_ref = set()

with gzip.open(strelka, 'rt') as strelka_data:
    for line in strelka_data:
        if line.startswith('#'):
            continue
        
        line = line.strip().split()
        
        chrom = line[0]
        pos = line[1]
        
        strelka_ref.add((chrom, pos))

with gzip.open(museq, 'rt') as museq_data, gzip.open(museq_filtered, 'wt') as museqout:
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