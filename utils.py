import errno
import gzip
import multiprocessing
import os
import pandas as pd

from subprocess import Popen, PIPE


class OpenFile(object):
    def __init__(self, filename, mode='rt'):
        self.filename = filename
        self.mode = mode

    def __enter__(self):
        if self.get_file_format(self.filename) in ["csv", 'plain-text']:
            self.handle = open(self.filename, self.mode)
        elif self.get_file_format(self.filename) == "gzip":
            self.handle = gzip.open(self.filename, self.mode)
        return self.handle

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.handle.close()

    def get_file_format(self, filepath):
        if filepath.endswith('.tmp'):
            filepath = filepath[:-4]

        _, ext = os.path.splitext(filepath)

        if ext == ".gz":
            return "gzip"
        else:
            return "plain-text"


def makedirs(directory, isfile=False):
    if isfile:
        directory = os.path.dirname(directory)
        if not directory:
            return

    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def build_shell_script(command, tag, tempdir):
    outfile = os.path.join(tempdir, "{}.sh".format(tag))
    with open(outfile, 'w') as scriptfile:
        scriptfile.write("#!/bin/bash\n")
        if isinstance(command, list) or isinstance(command, tuple):
            command = ' '.join(command) + '\n'
        scriptfile.write(command)
    return outfile


def run_cmd(cmd, output=None):
    stdout = PIPE
    if output:
        stdout = open(output, "w")

    cmd = ' '.join(cmd)
    print(cmd)

    p = Popen(cmd, stdout=stdout, stderr=PIPE, shell=True)

    cmdout, cmderr = p.communicate()
    retc = p.returncode

    if retc:
        raise Exception(
            "command failed. stderr:{}, stdout:{}".format(
                cmdout,
                cmderr))

    if output:
        stdout.close()


def run_in_gnu_parallel(commands, tempdir, ncores=None):
    makedirs(tempdir)

    scriptfiles = []

    for tag, command in enumerate(commands):
        scriptfiles.append(build_shell_script(command, tag, tempdir))

    parallel_outfile = os.path.join(tempdir, "commands")
    with open(parallel_outfile, 'w') as outfile:
        for scriptfile in scriptfiles:
            outfile.write("sh {}\n".format(scriptfile))

    if not ncores:
        ncores = str(multiprocessing.cpu_count())

    gnu_parallel_cmd = ['parallel', '--jobs', ncores, '<', parallel_outfile]

    run_cmd(gnu_parallel_cmd)


def filter_vcfs(sample_id, museq_vcf, strelka_vcf, work_dir):
    '''
    Original code by Diljot Grewal

    museq_paired and strekla_snv,
    take position intersection plus probability filter of 0.85
    (keep positions >= 0.85 that are in both)

    Modifications: take museq and strelka directly as input, output
    to temp_dir, take sample id as input and append to output
    filtered filename
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


def add_counts_to_maf(sample_id, temp_dir):
    n_counts = pd.read_csv(temp_dir + sample_id + '.csv', dtype=object, sep='\t').rename(columns={'chrom': 'Chromosome', 'coord': 'Start_Position', 'ref': 'Reference_Allele', 'alt': 'Tumor_Seq_Allele2'})
    
    maf = pd.read_csv(temp_dir + sample_id + '.maf', dtype=object, sep='\t', skiprows=1)
    maf = maf.drop(columns=['n_ref_count', 'n_alt_count'])
    maf = maf.merge(n_counts, on=['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2'], how='left')
    
    if os.path.exists(temp_dir + sample_id + '_tumour_counts.csv'):
        t_counts = pd.read_csv(temp_dir + sample_id + '_tumour_counts.csv', dtype=object, sep='\t').rename(columns={'chrom': 'Chromosome', 'coord': 'Start_Position', 'ref': 'Reference_Allele', 'alt': 'Tumor_Seq_Allele2'})
        maf = maf.drop(columns=['t_ref_count', 't_alt_count'])
        maf = maf.merge(t_counts, on=['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2'], how='left')
    
    maf.to_csv(temp_dir + sample_id + '.maf', index=None, sep='\t')
