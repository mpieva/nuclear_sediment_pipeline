import sys
from csv import reader
from Bio import SeqIO
from collections import defaultdict
import argparse
import os
#grep 'scientific name' ~/MetaGen/references/names.dmp |cut -d'|' -f 1,2>names_trimmed.dmp
#cut -d '|' -f 1,2,3 ~/MetaGen/references/naodes.dmp >nodes_trimmed.dmp

parser = argparse.ArgumentParser(description='Create a report from a kraken output.')
parser.add_argument('--db', required=True,
                    help='The kraken database to use')
parser.add_argument('--zeros', action='store_true',
                    help='Show also 0')
parser.add_argument('--clades', default=1, type=int,
                    help='Select only clade')
parser.add_argument('--minp', default=0.0, type=float,
                    help='Filter on the minimum percent of sequences for this clade')
parser.add_argument('--min', default=0, type=int,
                    help='Filter on the minimum sequences for this clade')
parser.add_argument('--rank',  help='Only return clades for specified rank')
parser.add_argument('--extractFile',  help='Fasta file where to extract seqence')
parser.add_argument('infile', metavar="kraken.output")

args = parser.parse_args()

db_prefix = os.path.abspath(args.db)
if args.rank and len(args.rank) > 1:
    args.rank = rank_code(rank)

seq_ids = defaultdict(list)
extract_ids = []

# remember to use filtered nodes.dmp and names.dmp
def load_taxonomy(db_prefix):
    name_map = {}
    rank_map = {}
    child_lists = defaultdict(list)
    with open(db_prefix+"/taxonomy/names_trimmed.dmp", 'r') as name_file:
        for line in name_file:
            node_id, name = line.strip().split('|')
            node_id = int(node_id)
            name = name.strip()
            name_map[node_id] = name
            
    with open(db_prefix+"/taxonomy/nodes_trimmed.dmp", 'r') as nodes_file:
        for line in nodes_file:
            node_id, parent_id, rank = line.strip().split('|')
            node_id = int(node_id)
            parent_id = int(parent_id)
            rank = rank.strip()
            if node_id == 1:
                parent_id = 0
            child_lists[parent_id].append(node_id)
            rank_map[node_id] = rank
    return (name_map, rank_map, child_lists)

def rank_code(rank):
    if rank == "species": return "S"
    if rank == "genus": return "G"
    if rank == "family": return "F"
    if rank == "order": return "O"
    if rank == "class": return "C"
    if rank == "phylum": return "P"
    if rank == "kingdom": return "K"
    if rank == "superkingdom": return "D"
    return "-"

def extract_seq_from_id(fileout, id_list, seqfile):
    num_seq_to_extract = len(id_list)
    with open(fileout+".fa", 'w') as fout:
        for rec in SeqIO.parse(seqfile, 'fasta'):
            if not num_seq_to_extract:
                break
            if rec.id in id_list:
                num_seq_to_extract -= 1
                #print(rec, file=sys.stderr)
                SeqIO.write(rec,fout,'fasta')

def dfs_report (node, depth):
    global extract_ids
    t_counts, c_counts, rank = taxo_counts[node], clade_counts[node], rank_map[node]
    if (not c_counts and not args.zeros):
        return
    c_counts_percent = round(c_counts * 100 / seq_count, 2)
    #filter on min seqences on clade
    #filter on min percent
    #filter on rank
    if not args.rank or (args.rank == rank_code(rank)) and (c_counts >= args.min and c_counts_percent >= args.minp):
        print ("{:6.2f}\t{}\t{}\t{}\t{}\t{}{}".format(
            c_counts_percent,
            c_counts,
            t_counts,
            rank_code(rank),
            node,
            "  " * depth,
            name_map[node]))
    # start saving the sequence mames for this clade
    if args.rank == rank_code(rank): extract_ids = []
    children = child_lists[node]
    if len(children):
        sorted_children = sorted(children, key=lambda k: clade_counts[k], reverse=True)
        #format output only if not filtered by rank
        if not args.rank : depth += 1
        for child in sorted_children:
            dfs_report(child, depth)
    if t_counts:# add only if the node has sequences assigned to it
        extract_ids.extend(seq_ids[node])
    # we want to extract up to a certain clade from a certain rank,
    # if there is a min sequences to extract, and only if a ref file is provided
    if (args.clades == node or rank_code(rank) == args.rank) and args.extractFile and len(extract_ids) >= args.min:
        print ("Extracting",len(extract_ids),"sequences for",name_map[node], file=sys.stderr)
        extract_seq_from_id(name_map[node], extract_ids, args.extractFile)
        extract_ids = []

def dfs_summation(node):
    children = child_lists[node]
    if len(children):
        for child in children:
            dfs_summation(child)
            clade_counts[node] += clade_counts.get(child, 0)

name_map, rank_map, child_lists = load_taxonomy(db_prefix)

print("dmp files loaded", file=sys.stderr)

seq_count = 0
taxo_counts = defaultdict(int)

with open(args.infile, 'r', newline='') as krakenfile:
    kfile = reader(krakenfile, delimiter='\t')
    for row in kfile:
        taxo_counts[int(row[2])] += 1
        seq_count += 1
        seq_ids[int(row[2])].append(row[1])

print(args.infile,"parsed", file=sys.stderr)

classified_count = seq_count - taxo_counts[0]
clade_counts = taxo_counts.copy()

#clade_counts

dfs_summation(args.clades)

unclassified_percent = 100
if seq_count:
    unclassified_percent = clade_counts.get(0) * 100 / seq_count

print ("{:6.2f}\t{}\t{}\t{}\t{}\t{}{}".format(
    unclassified_percent,
    clade_counts.get(0), taxo_counts[0], 
    "U", 0, "", "unclassified"))
dfs_report(args.clades, 0)
