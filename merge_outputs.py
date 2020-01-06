import glob
import pandas as pd

# merge multiple gistic or integer gene data text files
files_to_merge = glob.glob('merge/*.txt')
merged_file = open('merged.txt','w+')
first_file_to_merge = open(files_to_merge.pop(),'r')
for line in first_file_to_merge:
    merged_file.write(line)

merged_file.close()

merged_file = pd.read_csv('merged.txt', delimiter='\t', dtype={'entrez_id': str})
for file in files_to_merge:
    data_file = pd.read_csv(file, delimiter='\t', dtype={'entrez_id': str})
    merged_file = pd.merge(merged_file, data_file, on=['hugo_symbol', 'entrez_id'], how='outer')

merged_file.to_csv('merged.txt', index=None, sep='\t')
