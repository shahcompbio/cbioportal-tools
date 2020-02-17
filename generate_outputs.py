import click
import csv
import logging
import re

from io import StringIO
from scipy.stats import norm


# args mus be full filepaths
def extract(gtf, hgnc, igv_segs, titan_segs):
    # extract required columns from input files, place in single file
    extracted_file = StringIO()
    extracted_file.write('chr\tseg_start\tseg_end\tcopy_number\ttitan_state\tnum.mark\tmedian_logr\tensembl_id\thugo_symbol\tentrez_id\tgene_start\tgene_end\n')
    
    # test file: igv_segs.txt
    igv_file = open(igv_segs, 'r')
    next(igv_file)
    igv_reader = csv.reader(igv_file, delimiter='\t')
    
    # test file: titan_segs.txt
    titan_file = open(titan_segs, 'r')
    next(titan_file)
    titan_reader = csv.reader(titan_file, delimiter='\t')
    
    # test file: custom.txt
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

    # finally, write specified columns to extract.txt
    for igv_line, titan_line in zip(igv_reader, titan_reader):
        ensembl_ids = re.findall(r'ENSG\d+', titan_line[-1])
        # chr   seg_start   seg_end copy_number titan_state    num.mark    median_logr
        non_ensembl_info = titan_line[1] + '\t' + titan_line[2] + '\t' + titan_line[3] + '\t' + titan_line[9] + '\t' + titan_line[7] + '\t' + igv_line[4] + '\t' + titan_line[6]
        
        if not ensembl_ids:
            extracted_file.write(non_ensembl_info + '\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\t' + '' + '\n')
        else:
            for ensembl_id in ensembl_ids:
                extracted_file.write(non_ensembl_info + '\t' + ensembl_id + '\t')
                
                if ensembl_id in hgnc_dict:
                    # hugo_symbol   entrez_id
                    extracted_file.write(hgnc_dict[ensembl_id][0] + '\t' + hgnc_dict[ensembl_id][1] + '\t')
                else:
                    extracted_file.write('' + '\t' + '' + '\t')
                
                if ensembl_id in gtf_dict:
                    # gene_start    gene_end
                    extracted_file.write(gtf_dict[ensembl_id][0] + '\t' + gtf_dict[ensembl_id][1])
                else:
                    extracted_file.write('' + '\t' + '')

                extracted_file.write('\n')

    extracted_file.seek(0)
    return extracted_file 


def calculate_weighted_average(ensembl_dict, column_to_use):
    calculated_values = {}
    for ensembl_id in ensembl_dict:
        # find start and end points for all the segments gene is in
        seg_starts = [start for start in ensembl_dict[ensembl_id][0]]
        seg_ends = [end for end in ensembl_dict[ensembl_id][1]]
        values_to_use = [val for val in column_to_use[ensembl_id]]
        segs_to_remove = []
        
        # if ensembl_id gene is only present in one segment, add
        # associated copy number to calculated_values
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
        
        # determine remaining required information for calculation
        denominator_rest = 1 * len(seg_starts)
        numerator_rest = 0  
        for value in seg_starts:
            numerator_rest = numerator_rest + values_to_use[seg_starts.index(value)]

        # perform final weighted average calculation
        calculated_values[ensembl_id] = (numerator_start + numerator_rest + numerator_end) / (denominator_start + denominator_rest + denominator_end)

    return calculated_values


def transform(extracted_file, show_missing_hugo, show_missing_entrez, show_missing_both):
    '''
    perform weighted average calculations, and transformations
    
    ensembl_dict will store information for calculations
    gene_dict will store information for gene data output
    seg_dict will store information for segment data output
    '''
    
    ensembl_dict, gene_dict, seg_dict  = {}, {}, {}
    homd_segs, missing_hugo_symbol, missing_entrez_id, missing_both = [], [], [], []
    next(extracted_file)
    file_reader = csv.reader(extracted_file, delimiter='\t')
    for line in file_reader:    
        # if line has an associated ensembl_id and
        # segment doesn't have a length of zero
        if line[7] and line[1] != line[2]:
            # if gene is missing hugo_symbol or entrez_id
            if line[8] != '' or line[9] != '':
                if line[8] == '':
                    missing_hugo_symbol.append(line[7])
                if line[9] == '':
                    missing_entrez_id.append(line[7])
                
                # ensembl_id: [entrez_id, hugo_symbol]
                gene_dict[line[7]] = [line[9], line[8]]
                
                # set up a key-value pair where the key is
                # an ensembl id and the value is a list containing
                # the associated segment start points, end points,
                # copy numbers, titan states, gene start point, and
                # gene end point
                if line[7] not in ensembl_dict:
                    ensembl_dict[line[7]] = [[], [], [], [], 0, 0]
                
                ensembl_dict[line[7]][0].append(int(line[1]))
                ensembl_dict[line[7]][1].append(int(line[2]))
                ensembl_dict[line[7]][2].append(int(line[3]))
                ensembl_dict[line[7]][3].append(int(line[4]))
                ensembl_dict[line[7]][4] = int(line[10])
                ensembl_dict[line[7]][5] = int(line[11])

            # if gene is missing hugo_symbol and entrez_id
            if line[8] == '' and line[9] == '':
                missing_both.append(line[7])

        # (seg_start, seg_end): [chr, num.mark, titan_state, median_logr]
        seg_dict[(line[1], line[2])] = [line[0], line[5], line[4], line[6]]

    copy_numbers = {}
    for ensembl_id in ensembl_dict:
        # key: ensembl_id, value: associated copy number(s)
        copy_numbers[ensembl_id] = ensembl_dict[ensembl_id][2]

    calculated_cns = calculate_weighted_average(ensembl_dict, copy_numbers)
    
    titan_states = {}
    for ensembl_id in ensembl_dict:
        # key: ensembl_id, value: associated titan state(s)
        titan_states[ensembl_id] = ensembl_dict[ensembl_id][3]

    calculated_tss = calculate_weighted_average(ensembl_dict, titan_states)

    for ensembl_id in calculated_cns:
        # append weighted average of calculated copy number and
        # titan state, for each ensembl_id
        gene_dict[ensembl_id].append(str(calculated_cns[ensembl_id]))
        # remember to round to nearest integer for integer data
        gene_dict[ensembl_id].append(str(round(calculated_tss[ensembl_id])))

    # create a baseline (mu) for gene copy number transformation
    calc_cn_list = [calculated_cns[key] for key in calculated_cns]
    mu, _ = norm.fit(calc_cn_list)

    # perform required gene transformations on copy number
    for ensembl_id in calculated_cns:
        cn = calculated_cns[ensembl_id]
        
        if cn < 1:
            gene_dict[ensembl_id][2] = '-2'
            if cn < 0:
                logging.warning(f'{ensembl_id} has a calculated cn value lesser than 0')
        elif 1 <= cn <= mu-1:
            gene_dict[ensembl_id][2] = '-1'
        elif mu-1 < cn < mu+1:
            gene_dict[ensembl_id][2] = '0'
        elif mu+1 <= cn < 6:
            gene_dict[ensembl_id][2] = '1'
        elif cn >= 6:
            gene_dict[ensembl_id][2] = '2'

    if show_missing_hugo:
        print('Ensembl IDs missing HUGO symbols:')
        print(missing_hugo_symbol)
    if show_missing_entrez:
        print('Ensembl IDs missing Entrez IDs:')
        print(missing_entrez_id)
    if show_missing_both:
        print('Ensembl IDs missing both HUGO symbols and Entrez IDs:')
        print(missing_both, '\n')
    
    return gene_dict, seg_dict


def load(gene_dict, seg_dict, sample_id, output_dir, output_gistic_gene, output_integer_gene, output_log_seg, output_integer_seg):
    '''
    split generated file into four outputs

    gene_dict contains key-value pairs of
    ensembl_id:
    [entrez_id, hugo_symbol, transformed_calc_cn, calc_titan_state]
    
    seg_dict contains key-value pairs of
    (seg_start, seg_end):
    [chr, num.mark, titan_state, median_logr]
    '''
    
    gene_header = 'Hugo_Symbol\tEntrez_Gene_Id\t' + sample_id + '\n'
    segment_header = 'ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n'
    
    if output_gistic_gene:
        gistic_gene_data = open(output_dir + sample_id + '_gistic_gene_data.txt', 'w+')
        gistic_gene_data.write(gene_header)
        
        for ensembl_id in gene_dict:
            gistic_gene_data.write(gene_dict[ensembl_id][1] + '\t' + gene_dict[ensembl_id][0] + '\t' + gene_dict[ensembl_id][2] + '\n')
    
    if output_integer_gene:
        integer_gene_data = open(output_dir + sample_id + '_integer_gene_data.txt', 'w+')
        integer_gene_data.write(gene_header)

        for ensembl_id in gene_dict:
            integer_gene_data.write(gene_dict[ensembl_id][1] + '\t' + gene_dict[ensembl_id][0] + '\t' + gene_dict[ensembl_id][3] + '\n')
    
    if output_log_seg:
        log_seg_data = open(output_dir + sample_id + '_log_seg_data.txt', 'w+')
        log_seg_data.write(segment_header)

        for seg_length in seg_dict:
            log_seg_data.write(sample_id + '\t' + seg_dict[seg_length][0] + '\t' + seg_length[0] + '\t' + seg_length[1] + '\t' + seg_dict[seg_length][1] + '\t' + seg_dict[seg_length][3] + '\n')
    
    if output_integer_seg:
        integer_seg_data = open(output_dir + sample_id + '_integer_seg_data.txt', 'w+')
        integer_seg_data.write(segment_header)

        for seg_length in seg_dict:
            integer_seg_data.write(sample_id + '\t' + seg_dict[seg_length][0] + '\t' + seg_length[0] + '\t' + seg_length[1] + '\t' + seg_dict[seg_length][1] + '\t' + seg_dict[seg_length][2] + '\n')
    

@click.command()
@click.argument('gtf')
@click.argument('hgnc')
@click.argument('igv_segs')
@click.argument('titan_segs')
@click.argument('sample_id')
@click.option('--output_dir', default='')
@click.option('--output_gistic_gene/--no_gistic_gene', default=True)
@click.option('--output_integer_gene/--no_integer_gene', default=False)
@click.option('--output_log_seg/--no_log_seg', default=True)
@click.option('--output_integer_seg/--no_integer_seg', default=False)
@click.option('--show_missing_hugo/--no_missing_hugo', default=False)
@click.option('--show_missing_entrez/--no_missing_entrez', default=False)
@click.option('--show_missing_both/--no_missing_both', default=False)
def main(gtf, hgnc, igv_segs, titan_segs, sample_id, output_dir, output_gistic_gene, output_integer_gene, output_log_seg, output_integer_seg, show_missing_hugo, show_missing_entrez, show_missing_both):
    extracted_file = extract(gtf, hgnc, igv_segs, titan_segs)
    gene_dict, seg_dict = transform(extracted_file, show_missing_hugo, show_missing_entrez, show_missing_both)
    load(gene_dict, seg_dict, sample_id, output_dir, output_gistic_gene, output_integer_gene, output_log_seg, output_integer_seg)


if __name__ == '__main__':
    main()
