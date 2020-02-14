## cBioPortal Copy Number Data Formatting Tools

### generate_outputs.py

Code that takes four files as input:

1. genenames.org custom text file download with Approved symbol (HUGO), NCBI Gene ID (Entrez Gene ID), and Ensembl gene ID: `custom.txt`
2. a .gtf text file containing human gene information (gene start, end, and id information is used): `Homo_sapiens.GRCh37.73.gtf` (not available in the repo due to GitHub file size constraints)
3. Integrative Genomics Viewer segments text file: `igv_segs.txt`
4. TITAN segments text file: `titan_segs.txt`

and returns four text files:

1. gistic format segment data
2. gistic format gene data
3. integer segment data
4. integer gene data

For testing please run `python test_generate.py test/generate_outputs/test_input/Homo_sapiens.GRCh37.73.gtf test/generate_outputs/test_input/custom.txt test/generate_outputs/test_input/igv_segs.txt test/generate_outputs/test_input/titan_segs.txt test_sample_id`

### merge_outputs.py

Code that takes a directory of gistic format gene data OR integer gene data text files and outputs a single text file with their contents merged.

For testing please run `python test_merge.py`
