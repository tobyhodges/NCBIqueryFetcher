#! /usr/bin/python

import xml.etree.ElementTree as ET
import urllib2
import getpass
from optparse import OptionParser
from time import sleep
from Bio import SeqIO

##############################################################################
# FUNCTIONS                                                                  #
##############################################################################

# set proxy settings to allow entrez queries to be sent through a proxy server
# using Eutils
def buildProxySettings(serverAddress):
  usrnm = getpass.getuser() # get username from session login details
  pswd = getpass.getpass('please enter your password for the proxy server: ') # ask user for their proxy server password
  proxy = urllib2.ProxyHandler({'http': ''.join(['http://', usrnm,':',pswd,'@', serverAddress])}) # build proxy server settings for urllib2
  opener = urllib2.build_opener(proxy) # create proxy-aware urllib2 opener
  return opener

# get the list of IDs from a provided file. IDs must be one-per-line
# in the file, with no empty lines
def getRecordIDsFromFile(IDlistFilename):
  IDlist = []
  IDFH = open(IDlistFilename, 'r')
  for line in IDFH.readlines():
    line = line.strip()
    IDlist.append(line)
  return IDlist

# takes a database name and a dictionary of search fields and terms to be used in the search
# returns a list of record IDs that can be used in a second Eutils query
def getRecordIDsFromNCBI(database, searchTermDictionary):
  baseURL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/' # set Eutils base URL
  extension = 'esearch.fcgi?db=' # start building search term extension
  extension += database # specify database to be searched
  extension += '&retmax=100000'
  extension += '&term=' # get ready to add terms
  termCount = 0
  for searchTerm in searchTermDictionary: # add search terms to extension
    if termCount > 0:
      extension += '+AND+'
    termValue = searchTermDictionary[searchTerm]
    termValue = termValue.replace(' ', '+')
    termString = termValue + '[' + searchTerm.replace(' ', '+') + ']'
    extension += termString
  getIDsURL = baseURL + extension # build full search URL
  IDconnx = urllib2.urlopen(getIDsURL) # open conenction to NCBI
  IDtree = ET.parse(IDconnx) # get search output
  IDroot = IDtree.getroot()
  IDlist = []
  for child in IDroot: # build list of IDs returned from search
    if child.tag == 'IdList':
      for Id in child.findall('Id'):
        IDlist.append(Id.text)
  return IDlist

# takes a list of nucleotide database ID strings and fetches the associated sequence records
# returns the sequences as a list of Biopython Sequence Records
def getRecords(IDlist, format):
  if len(IDlist) < 500:
    baseURL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    extension = 'efetch.fcgi?db=nuccore&id='
    IDsAsString = ','.join(IDlist)
    extension += IDsAsString
    extension += '&rettype='
    extension += format
    extension += '&retmode=text'
    getRecordsURL = baseURL + extension
    RecordsConnx = urllib2.urlopen(getRecordsURL)
    seqRecords = SeqIO.parse(RecordsConnx, format)
  else:
    seqRecords = []
    chunkSize = 500
    for i in xrange(0, len(IDlist), chunkSize):
      baseURL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
      extension = 'efetch.fcgi?db=nuccore&id='
      IDblock = IDlist[i:i+chunkSize]
      IDsAsString = ','.join(IDblock)
      extension += IDsAsString
      extension += '&rettype='
      extension += format
      extension += '&retmode=text'
      getRecordsURL = baseURL + extension
      RecordsConnx = urllib2.urlopen(getRecordsURL)
      seqRecordBlock = list(SeqIO.parse(RecordsConnx, format))
      seqRecords += seqRecordBlock
      sleep(5)
  return seqRecords

##############################################################################
#                                                                            #
##############################################################################
# SCRIPT OPTIONS                                                             #
##############################################################################
# set script options, allowing the user to specify which fields they would
# like to search against, and the values associated with these fields
searchTermParser = OptionParser(usage = "%prog [options] <output_file.fasta>", \
                      description = "fetches FASTA format records from NCBI nucleotide database.")
searchTermParser.add_option('-a', '--accession', dest='accession', type='string', help='nucleotide database accession', default='', metavar='STR')
searchTermParser.add_option('-x', '--all', dest='all', type = 'string', help='all fields', default='', metavar='STR')
searchTermParser.add_option('-u', '--author', dest='author', type='string', help='author associated with submission of record', default='', metavar='STR')
searchTermParser.add_option('-e', '--ec_number', dest='ecnum', type='string', help='enzyme commission number associated with the record', default = '', metavar = 'STR')
searchTermParser.add_option('-k', '--feature_key', dest='fkey', type='string', help='key in feature table', default = '', metavar = 'STR')
searchTermParser.add_option('-f', '--filter', dest='filter', type='string', help='a filtered subset of the database', default = '', metavar = 'STR')
searchTermParser.add_option('-g', '--gene_name', dest='gene', type='string', help='gene name associated with record', default = '', metavar = 'STR')
searchTermParser.add_option('-r', '--genome_project', dest='proj', type='string', help='numeric id of genome project', default = '', metavar = 'STR')
searchTermParser.add_option('-j', '--journal', dest='journal', type='string', help='journal in which the record was originally published', default = '', metavar = 'STR')
searchTermParser.add_option('-w', '--keyword', dest='keyword', type='string', help='keyword in entry', default = '', metavar = 'STR')
searchTermParser.add_option('-m', '--modification_date', dest='moddate', type='string', help='date of most recent modification to record', default = '', metavar = 'STR')
searchTermParser.add_option('-o', '--organism', dest='organism', type='string', help='organism with which the record is associated', default = '', metavar = 'STR')
searchTermParser.add_option('-q', '--primary_accession', dest='primacc', type='string', help='first accession associated with record', default = '', metavar = 'STR')
searchTermParser.add_option('-c', '--properties', dest='prop', type='string', help='properties associated with record (molecular type, source database, etc.)', default = '', metavar = 'STR')
searchTermParser.add_option('-p', '--protein_name', dest='protein', type='string', help='name of protein associated with record', default = '', metavar = 'STR')
searchTermParser.add_option('-d', '--publication_date', dest='date', type='string', help='date that record was first published on Entrez', default = '', metavar = 'STR')
searchTermParser.add_option('-i', '--sequence_id', dest='seqid', type='string', help='NCBI identifier', default = '', metavar = 'STR')
searchTermParser.add_option('-l', '--sequence_length', dest='length', type='string', help='total length of sequence (can be given as a range in format \'a:b\')', default = '', metavar = 'STR')
searchTermParser.add_option('-n', '--substance_name', dest='substance', type='string', help='chemical substances associated with record', default = '', metavar = 'STR')
searchTermParser.add_option('-t', '--title', dest='title', type='string', help='in title of record', default = '', metavar = 'STR')
searchTermParser.add_option('--format', dest='outputFormat', type='string', help='desired format of output file (default=fasta, for full genbank format use \'gb\')', default='fasta', metavar='STR')
searchTermParser.add_option('--id_file', dest='idListFile', type='string', help='a list of IDs to obtain records for, if you wish to forego the search process', default='', metavar='STR')
searchTermParser.add_option('--proxy', dest='proxyServer', type='string', help='send queries through proxy server at this address. Specify this option if you are working behind a proxy server.', default=None, metavar='STR')

##############################################################################
#
##############################################################################

if __name__ == '__main__':
  # parse options and argument
  (options, arguments) = searchTermParser.parse_args()
# assign output filename from parsed argument
  outputFile = arguments[0]
# get output file format from option field
  outputFormat = options.outputFormat
  if options.proxyServer != None:
# get proxy settings
    proxyOpener = buildProxySettings(options.proxyServer)
    urllib2.install_opener(proxyOpener)
# if a list of IDs has been provided, extract the IDs from the file
  if options.idListFile != '':
    idListFile = options.idListFile
    IDs = getRecordIDsFromFile(idListFile)
    print IDs
# otherwise, get the search terms from the option fields
  else:
    allSearchTerms = {}
    allSearchTerms['ACCN'] = options.accession
    allSearchTerms['ALL'] = options.all
    allSearchTerms['AUTH'] = options.author
    allSearchTerms['ECNO'] = options.ecnum
    allSearchTerms['FKEY'] = options.fkey
    allSearchTerms['FILT'] = options.filter
    allSearchTerms['GENE'] = options.gene
    allSearchTerms['Genome Project'] = options.proj
    allSearchTerms['JOUR'] = options.journal
    allSearchTerms['KYWD'] = options.keyword
    allSearchTerms['MDAT'] = options.moddate
    allSearchTerms['ORGN'] = options.organism
    allSearchTerms['PACC'] = options.primacc
    allSearchTerms['PROP'] = options.prop
    allSearchTerms['PROT'] = options.protein
    allSearchTerms['PDAT'] = options.date
    allSearchTerms['SQID'] = options.seqid
    allSearchTerms['SLEN'] = options.length
    allSearchTerms['SUBS'] = options.substance
    allSearchTerms['TITL'] = options.title
# filter the full set of search fields and include any non-empty fields
# in the set of terms to be searched
    searchTerms = {}
    for searchTerm in allSearchTerms:
      if allSearchTerms[searchTerm] != '':
        searchTerms[searchTerm] = allSearchTerms[searchTerm]
# check that user has entered at least one field to be searched
    if len(searchTerms.keys()) == 0:
      die('if you are not providing a list of IDs yourself, you must enter at least one search term!')
    else:
# fetch record IDs for search results
      IDs = getRecordIDsFromNCBI('nucleotide', searchTerms)
# get FASTA records associated with search results
  sequences = getRecords(IDs, outputFormat)
# write sequence records to output file
  OUTFH = open(outputFile, 'w')
  SeqIO.write(sequences, OUTFH, outputFormat)
