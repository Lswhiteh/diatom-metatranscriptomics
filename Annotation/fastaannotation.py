#!/usr/bin/env python3

from Bio import SeqIO
import csv
import sys

'''
Annotates Trinity fasta outputs with phylodb hits using DIAMOND tabular output

Usage:
    python3 fastqannotation.py <trinitycontigs.fasta> <diamondoutput.m8> <annotatedout.fasta> <mappingfile.tsv>
'''

# PhyloDB paths can be soft-coded if you comment out the file paths and uncomment the next two lines
# If this is done, the taxonomy and gene file will need to be supplied to the program call

taxonomy_file = "/pine/scr/l/s/lswhiteh/phylodbannotation/phylodb_1.076.taxonomy.txt"
gene_file = "/pine/scr/l/s/lswhiteh/phylodbannotation/phylodb_1.076.annotations.txt"
#taxonomy_file = sys.argv[5]
#annotated_file = sys.argv[6]

# Get data files from command line
contigs_file = sys.argv[1]
diamond_file = sys.argv[2]
annotated_file = sys.argv[3]
mapping_file = sys.argv[4]

# Load phylodb master-lists as dictionaries for easy parsing
diamond_dict = {}
tax_dict = {}
gene_dict = {}

with open(diamond_file) as diamondfile:
    for line in diamondfile:
        row = line.strip().split("\t")
        diamond_dict[row[0]] = row[1:]

with open(taxonomy_file) as taxfile:
    for line in taxfile:
        row = line.strip().split("\t")
        tax_dict[row[0]] = row[1:]

with open(gene_file) as genefile:
    for line in genefile:
        row = line.strip().split("\t")
        gene_dict[row[0]] = row[1:]

'''
1. Get each fasta descriptor
2. If descriptor is in diamond output, get diamond annotation
3. Use diamond annotation to get gene annotation from gene_dict
4. Use gene annotation as key for tax_dict
5. Write all fields in fasta descriptor as new header
'''




with open(contigs_file, 'rU') as contigs, open(annotated_file, 'w') as annotated, open(mapping_file, 'w') as mapping:
    records = SeqIO.parse(contigs, "fasta")
    writer = csv.writer(mapping, delimiter='\t')
    writer.writerow(["TrinityID", "Gene", "Organism", "Taxonomy"])
    for record in records:  
        if record.id in diamond_dict:
            
            old_record = record.id

            gene = diamond_dict[record.id][0]

            organism = gene_dict[gene][1]

            tax = tax_dict[organism][1]

            record.id = "\t".join([old_record, str(gene), organism, str(tax)])
            SeqIO.write(record, annotated, 'fasta-2line')

            writer.writerow([old_record, gene, organism, tax])