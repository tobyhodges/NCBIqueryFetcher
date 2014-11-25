# NCBIqueryFetcher.py
================

## Description

Fetch all entries in NCBI nucleotide database that match search criteria. 

This script uses the NCBI EUtils to compile a list of all nucleotide entries that match a user-defined set of search criteria, and then fetches these entries in the desired format.

The script splits fetching of records for large lists of GIs into smaller chunks, with a pause built in between fetches. This is designed to prevent the user from breaching the access limits defined by NCBI [here](http://www.ncbi.nlm.nih.gov/books/NBK25497/#_chapter2_Usage_Guidelines_and_Requiremen_). If you edit these settings, please make sure that you don't breach the access limits described at the previous link.

Also take note of the [NCBI disclaimer](www.ncbi.nlm.nih.gov/About/disclaimer.html) associated with use of the EUtils.

================

## Usage

**Usage: python NCBIqueryFetcher.py [options] <output_file.fasta>**  

fetches FASTA format records from NCBI nucleotide database.  

Options:  
  -h, --help                            show this help message and exit  
  -a STR, --accession=STR               nucleotide database accession  
  -x STR, --all=STR                     all fields  
  -u STR, --author=STR                  author associated with submission of record  
  -e STR, --ec_number=STR               enzyme commission number associated with the record  
  -k STR, --feature_key=STR             key in feature table  
  -f STR, --filter=STR                  a filtered subset of the database  
  -g STR, --gene_name=STR               gene name associated with record  
  -r STR, --genome_project=STR          numeric id of genome project  
  -j STR, --journal=STR                 journal in which the record was originally published  
  -w STR, --keyword=STR                 keyword in entry  
  -m STR, --modification_date=STR       date of most recent modification to record  
  -o STR, --organism=STR                organism with which the record is associated  
  -q STR, --primary_accession=STR       first accession associated with record  
  -c STR, --properties=STR              properties associated with record (molecular type, source database, etc.)  
  -p STR, --protein_name=STR            name of protein associated with record  
  -d STR, --publication_date=STR        date that record was first published on Entrez  
  -i STR, --sequence_id=STR             NCBI identifier  
  -l STR, --sequence_length=STR         total length of sequence (can be given as a range in format 'a:b')  
  -n STR, --substance_name=STR          chemical substances associated with record  
  -t STR, --title=STR                   in title of record  
  --format=STR                          desired format of output file (default=fasta, for full genbank format use 'gb')  
  --id_file=STR                         a list of IDs to obtain records for, if you wish to forego the search process  
  --proxy=STR                           send queries through proxy server at this address. Specify this option if you are working behind a proxy server.
