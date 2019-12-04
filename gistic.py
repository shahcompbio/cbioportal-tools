import click
import csv
import re


# files need to live in the same directory as script, arg mus be full filename
def extract(gtf, hgnc, igv_segs, titan_segs):
	# extract required information from input files, place in single file
	extract_file = open("extract.txt","w+")
	extract_file.write("chr\tseg_start\tseg_end\tmedian_logr\ttitan_state\tnum.mark\tensembl_id\thugo_symbol\tentrez_id\tgene_start\tgene_end\n")
	
	# from igv_segs.txt:
	# num.mark
	igv_file = open(igv_segs, 'r')
	next(igv_file)
	igv_reader = csv.reader(igv_file, delimiter='\t')
	
	# from titan_segs.txt:
	# Chromosome, Start_Position(bp), End_Position(bp), Median_logR,
	# TITAN_state, Pygenes(gene_id,gene_name;) (throw away gene_name)
	titan_file = open(titan_segs, 'r')
	next(titan_file)
	titan_reader = csv.reader(titan_file, delimiter='\t')
	
	# from custom.txt
	# Approved symbol, NCBI Gene ID
	hgnc_file = open(hgnc, 'r')
	next(hgnc_file)
	hgnc_reader = csv.reader(hgnc_file, delimiter='\t')
	hgnc_dict = {}
	for hgnc_line in hgnc_reader:
		if len(hgnc_line) == 3:
			hgnc_dict[hgnc_line[2]] = (hgnc_line[0], hgnc_line[1])
	
	# use Homo_sapiens.GRCh37.73.gtf.txt columns 4 and 5
	# (gene start and end), for associated gene_id
	gtf_file = open(gtf, 'r')
	gtf_reader = csv.reader(gtf_file, delimiter='\t')
	gtf_dict = {}
	# TODO: extract relevant gene start and end information

	for igv_line, titan_line in zip(igv_reader, titan_reader):
		ensembl_ids = re.findall(r'ENSG\d+', titan_line[-1])
		non_ensembl_info = titan_line[1] + '\t' + titan_line[2] + '\t' + titan_line[3] + '\t' + titan_line[6] + '\t' + titan_line[7] + '\t' + igv_line[4]
		if not ensembl_ids:
			extract_file.write(non_ensembl_info + '\n')
		else:
			for ensembl_id in ensembl_ids:
				extract_file.write(non_ensembl_info + '\t' + ensembl_id)
				if ensembl_id in hgnc_dict:
					extract_file.write('\t' + hgnc_dict[ensembl_id][0] + '\t' + hgnc_dict[ensembl_id][1])
				# if ensembl_id in gtf_dict:
				# 	extract_file.write('\t' + gtf_dict[ensembl_id][0] + '\t' + gtf_dict[ensembl_id][1])

				extract_file.write('\n')

	extract_file.close()
	return extract_file


def transform():
	# perform weighted averge calculations and required transformations
	pass


def load():
	# split generated file into the four output 
	pass


@click.command()
@click.argument('gtf')
@click.argument('hgnc')
@click.argument('igv_segs')
@click.argument('titan_segs')
@click.argument('sample_id')
def main(gtf, hgnc, igv_segs, titan_segs, sample_id):
    files = extract(gtf, hgnc, igv_segs, titan_segs)
    transform()
    load()


if __name__ == '__main__':
    main()