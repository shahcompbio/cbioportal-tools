## cBioPortal Copy Number Formatting
This repository contains code that takes four files as input: 
1. genenames.org custom text file download with Approved symbol (HUGO), NCBI Gene ID (Entrez Gene ID), and Ensembl gene ID: `custom.txt`
2. a .gtf text file containing human gene information (gene start, end, and id information is used): `Homo_sapiens.GRCh37.73.gtf`
3. Integrative Genomics Viewer segments text file: `igv_segs.txt`
4. TITAN segments text file: `titan_segs.txt` 

and returns four text files:  
1. gistic format segment data
2. gistic format gene data
3. integer segment data
4. integer gene data
