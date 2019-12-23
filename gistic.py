import click
import csv
import re

from scipy.stats import norm


# files need to live in the same directory as script, arg mus be full filename
def extract(gtf, hgnc, igv_segs, titan_segs):
	# extract required information from input files, place in single file
	extracted_file = open('extract.txt','w+')
	extracted_file.write('chr\tseg_start\tseg_end\tcopy_number\ttitan_state\ttitan_call\tnum.mark\tensembl_id\thugo_symbol\tentrez_id\tgene_start\tgene_end\n')
	
	# from igv_segs.txt:
	# num.mark
	igv_file = open(igv_segs, 'r')
	next(igv_file)
	igv_reader = csv.reader(igv_file, delimiter='\t')
	
	# from titan_segs.txt:
	# Chromosome, Start_Position(bp), End_Position(bp), Copy_Number,
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
	ensembl_set = set()
	for gtf_line in gtf_reader:
		ensembl_id = re.findall(r'ENSG\d+', gtf_line[-1])[0]
		ensembl_set.add(ensembl_id)

	gtf_file.seek(0)
	for ensembl_id in ensembl_set:
		gtf_dict[ensembl_id] = [[], []]

	# find all start and end points for each particular gene
	for gtf_line in gtf_reader:
		ensembl_id = re.findall(r'ENSG\d+', gtf_line[-1])[0]
		gtf_dict[ensembl_id][0].append(gtf_line[3])
		gtf_dict[ensembl_id][1].append(gtf_line[4])

	# consolidate sections by taking one start and one end coordinate
	for ensembl_id in gtf_dict:
		gtf_dict[ensembl_id][0] = min(gtf_dict[ensembl_id][0])
		gtf_dict[ensembl_id][1] = max(gtf_dict[ensembl_id][1])

	for igv_line, titan_line in zip(igv_reader, titan_reader):
		ensembl_ids = re.findall(r'ENSG\d+', titan_line[-1])
		non_ensembl_info = titan_line[1] + '\t' + titan_line[2] + '\t' + titan_line[3] + '\t' + titan_line[9] + '\t' + titan_line[7] + '\t' + titan_line[8] + '\t' + igv_line[4]
		if not ensembl_ids:
			extracted_file.write(non_ensembl_info + '\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\n')
		else:
			for ensembl_id in ensembl_ids:
				extracted_file.write(non_ensembl_info + '\t' + ensembl_id + '\t')
				if ensembl_id in hgnc_dict:
					extracted_file.write(hgnc_dict[ensembl_id][0] + '\t' + hgnc_dict[ensembl_id][1] + '\t')
				else:
					extracted_file.write('' + '\t' + '' + '\t')
				
				if ensembl_id in gtf_dict:
					extracted_file.write(gtf_dict[ensembl_id][0] + '\t' + gtf_dict[ensembl_id][1])
				else:
					extracted_file.write('' + '\t' + '')

				extracted_file.write('\n')

	extracted_file.close() 


def calculate_weighted_average(ensembl_dict, column_to_use):
	calculated_values = {}
	for ensembl_id in ensembl_dict:
		# find start and end points for all the segments gene is in
		seg_starts = [start for start in ensembl_dict[ensembl_id][0]]
		seg_ends = [end for end in ensembl_dict[ensembl_id][1]]
		values_to_use = [val for val in column_to_use[ensembl_id]]
		segs_to_remove = []
		
		# check if segment has a length of 0
		for seg_start, seg_end in zip(seg_starts, seg_ends):
			if seg_start == seg_end:
				segs_to_remove.append(seg_start)
		
		# if it does, remove it
		for seg in segs_to_remove:
			del values_to_use[seg_starts.index(seg)]
			seg_starts.remove(seg), seg_ends.remove(seg)
		
		# if ensembl_id gene is only present in one segment, add
		# associated copy number to calculated_values. If it is in a 0-length 
		# segment, and thus seg_starts and seg_ends are empty, move
		# to next ensembl_id
		if len(seg_starts) <= 1:
			if len(seg_starts) == 1:
				calculated_values[ensembl_id] = values_to_use[0]
			
			continue

		gene_start = ensembl_dict[ensembl_id][4]
		gene_end = ensembl_dict[ensembl_id][5]
		
		denominator_start = (min(seg_ends) - gene_start) / (min(seg_ends) - min(seg_starts))
		denominator_end = ((max(seg_ends) - max(seg_starts)) - (max(seg_ends) - gene_end)) / (max(seg_ends) - max(seg_starts))
		numerator_start = denominator_start * values_to_use[seg_starts.index(min(seg_starts))] 
		numerator_end = denominator_end * values_to_use[seg_starts.index(max(seg_starts))]
		
		# remove calculation values from used segs
		values_to_remove = [values_to_use[seg_starts.index(min(seg_starts))], values_to_use[seg_starts.index(max(seg_starts))]]
		for value in values_to_remove:
			values_to_use.remove(value)
		
		# remove min and max segment start and end coordinates
		# this is to easily iterate over remaining segments
		seg_starts.remove(min(seg_starts)), seg_starts.remove(max(seg_starts))
		seg_ends.remove(min(seg_ends)), seg_ends.remove(max(seg_ends))
		
		# calculate remaining required information
		denominator_rest = 1 * len(seg_starts)
		numerator_rest = 0	
		for value in seg_starts:
			numerator_rest = numerator_rest + values_to_use[seg_starts.index(value)]

		# perform final weighted average calculation
		calculated_values[ensembl_id] = (numerator_start + numerator_rest + numerator_end) / (denominator_start + denominator_rest + denominator_end)

	return calculated_values


def transform():
	# perform weighted average calculations, and transformations
	# ensembl_dict will store information for calculations
	# gene_dict will store information for gene data output
	# seg_dict will store information for segment data output
	ensembl_dict, gene_dict, seg_dict  = {}, {}, {}
	homd_segs = []
	extracted_file = open('extract.txt','r')
	next(extracted_file)
	file_reader = csv.reader(extracted_file, delimiter='\t')
	for line in file_reader:
		# if line has an associated ensembl_id
		if line[7]:
			# ensembl_id = [entrez_id, hugo_symbol]
			gene_dict[line[7]] = [line[8], line[9]]

			# set up a key-value pair where the key is
			# an ensembl id and the value is a list containing
			# the associated segment start points, end points,
			# copy numbers, titans_states gene start point, and
			# gene end point
			if line[7] not in ensembl_dict:
				ensembl_dict[line[7]] = [[], [], [], [], 0, 0]
			
			ensembl_dict[line[7]][0].append(int(line[1]))
			ensembl_dict[line[7]][1].append(int(line[2]))
			ensembl_dict[line[7]][2].append(float(line[3]))
			ensembl_dict[line[7]][3].append(float(line[4]))
			ensembl_dict[line[7]][4] = int(line[10])
			ensembl_dict[line[7]][5] = int(line[11])

		# (seg_start, seg_end) = [chr, num.mark, titan_state, copy_number]
		seg_dict[(line[1], line[2])] = [line[0], int(line[6]), int(line[4]), int(line[3])]
		
		# check for homozygous deletion
		if line[5] == 'HOMD':
			homd_segs.append((line[1], line[2]))

	copy_numbers = {}
	for ensembl_id in ensembl_dict:
		copy_numbers[ensembl_id] = ensembl_dict[ensembl_id][2]

	calculated_cns = calculate_weighted_average(ensembl_dict, copy_numbers)
	
	titan_states = {}
	for ensembl_id in ensembl_dict:
		titan_states[ensembl_id] = ensembl_dict[ensembl_id][3]

	calculated_tss = calculate_weighted_average(ensembl_dict, titan_states)

	for ensembl_id in calculated_cns:
		gene_dict[ensembl_id].append(calculated_cns[ensembl_id])
		gene_dict[ensembl_id].append(calculated_tss[ensembl_id])

	# create a baseline (mu) for gene copy number transformation
	calc_cn_list = [calculated_cns[key] for key in calculated_cns]
	mu, std = norm.fit(calc_cn_list)

	# perform required gene transformations on copy number
	for ensembl_id in calculated_cns:
		cn = calculated_cns[ensembl_id]
		if 0 <= cn < 1:
			gene_dict[ensembl_id][2] = -2
		elif 1 <= cn <= mu-1:
			gene_dict[ensembl_id][2] = -1
		elif mu-1 < cn < mu+1:
			gene_dict[ensembl_id][2] = 0
		elif mu+1 <= cn < 6:
			gene_dict[ensembl_id][2] = 1
		elif cn >= 6:
			gene_dict[ensembl_id][2] = 2

	# create a baseline (mu) for segment copy number transformation
	cn_list = [seg_dict[ensembl_id][3] for ensembl_id in seg_dict]
	mu, std = norm.fit(cn_list)

	# perform required segment transformations on copy number
	for seg_length in seg_dict:
		cn = seg_dict[seg_length][3]
		if seg_length in homd_segs:
			seg_dict[seg_length][3] = -2
		elif 0 <= cn <= mu-1:
			seg_dict[seg_length][3] = -1
		elif mu-1 < cn < mu+1:
			seg_dict[seg_length][3] = 0
		elif mu+1 <= cn < 6:
			seg_dict[seg_length][3] = 1
		elif cn >= 6:
			seg_dict[seg_length][3] = 2

	return gene_dict, seg_dict


def load(gene_dict, seg_dict, sample_id):
	# split generated file into the four outputs
	gene_header = 'entrez_id\thugo_symbol\t' + sample_id + '\n'
	segment_header = 'sample_id\tchr\tseg_start\tseg_end\tnum.mark\tseg.mean\n'
	extracted_file = open('extract.txt','r')
	next(extracted_file)
	gistic_gene_data = open('gistic_gene_data.txt','w+')
	gistic_gene_data.write(gene_header)
	

@click.command()
@click.argument('gtf')
@click.argument('hgnc')
@click.argument('igv_segs')
@click.argument('titan_segs')
@click.argument('sample_id')
def main(gtf, hgnc, igv_segs, titan_segs, sample_id):
    extract(gtf, hgnc, igv_segs, titan_segs)
    gene_dict, seg_dict = transform()
    load(gene_dict, seg_dict, sample_id)


if __name__ == '__main__':
    main()