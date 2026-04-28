#!/usr/bin/env python
import re
import os, argparse

# sample_id = "ERR552113"
# metrics_file = "type_data/clean/ERR552113_readMetrics.tsv"
# quast_file = "type_data/quast/ERR552113_report.txt"
# kraken_file = "type_data/kraken/ERR552113_top_kraken_species_results"

parser = argparse.ArgumentParser()
parser.add_argument("-m", help="path to read metrics file")
parser.add_argument("-q", help="path to Quast file")
parser.add_argument("-k", help="path to Kraken file")
parser.add_argument("-i", help="the sample id")
# parser.add_argument("-u", help="username of the submitter")
args = parser.parse_args()

if not args.m or not os.path.isfile(args.m):
    print("Please specify a valid read metrics file path using the '-m' parameter")
    exit()
    
if not args.q or not os.path.isfile(args.q):
    print("Please specify a valid Quast file path using the '-q' parameter")
    exit()
    
if not args.k or not os.path.isfile(args.k):
    print("Please specify a valid Kraken metrics file path using the '-k' parameter")
    exit()
    
if not args.i:
    print("Please specify a valid sample id using the '-i' parameter")
    exit()
    
info_str = ""
# parse the readMetrics file; the second line contains all the desired fields
with open(args.m) as file:
    i = 0
    for metrics_str in file:
        metrics_str = metrics_str.replace("\n", "")
        if i == 0:
            i = i + 1
        elif i == 1:
            break

contigs = "N/A"
largest_contig = "N/A"
total_length = "N/A"
N50 = "N/A"
L50 = "N/A"
# parse the quast report file; there are separate lines for various fields e.g. # contigs; here "\s+" 
# is used as a separator because there can be many spaces between the name of the field and its value;
# also sometimes [2] is used with re.split instead of [1] because some fields have a space in their name
# e.g. "Largest contig"
with open(args.q) as file:
    for quast_str in file:
        if quast_str.startswith("# contigs") and "(" not in metrics_str:
            contigs = re.split("\s+", quast_str)[2]
        elif quast_str.startswith("Largest contig"):
            largest_contig = re.split("\s+", quast_str)[2]
        elif quast_str.startswith("Total length"):
            total_length = re.split("\s+", quast_str)[2]
        elif quast_str.startswith("N50"):
            N50 = re.split("\s+", quast_str)[1]
        elif quast_str.startswith("L50"):
            L50 = re.split("\s+", quast_str)[1]

species = "N/A"
# parse the kraken file; the second line contains the species name
with open(args.k) as file:
    i = 0
    for species_str in file:
        if i == 0:
            i = i + 1
        elif i == 1:
            species_str = species_str.replace("\n", "")
            species = re.split("\t", species_str)[5]
            species = species.strip()
            break

info_str = metrics_str + "\t" + contigs + "\t" + largest_contig + "\t" + total_length + "\t" + N50 + "\t" + L50 + "\t" +species + "\n"

text_file = open(args.i+"_isolate_info_file.tsv", "w")

#write string to file
text_file.write(info_str)

#close file
text_file.close()