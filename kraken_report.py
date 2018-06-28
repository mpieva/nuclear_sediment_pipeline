import sys
from csv import reader
from Bio import SeqIO
from pysam import AlignmentFile
from collections import defaultdict
import argparse
import os
import json
#grep 'scientific name' ~/MetaGen/references/names.dmp |cut -d'|' -f 1,2>names_trimmed.dmp
#cut -d '|' -f 1,2,3 ~/MetaGen/references/naodes.dmp >nodes_trimmed.dmp

parser = argparse.ArgumentParser(description='Create a report from a kraken output. Optionally extract reads')
parser.add_argument('--db', required=True,
                    help='The kraken database to use')
parser.add_argument('--zeros', action='store_true',
                    help='Show also 0')
parser.add_argument('--clades', default=False,
                    help='Select only specified clades (comma separated)')
parser.add_argument('--minp', default=0.0, type=float,
                    help='Filter on the minimum percent of sequences for this clade')
parser.add_argument('--min', default=0, type=int,
                    help='Filter on the minimum sequences for this clade')
parser.add_argument('--rank',  help='Only return clades for specified rank')
parser.add_argument('--extractFile',  help='File where to extract sequence from')
parser.add_argument('infile', metavar="kraken.output")

args = parser.parse_args()

if args.clades: # handle providing multiple clades, comma separated
    args.clades = args.clades.split(",")

db_prefix = os.path.abspath(args.db)
if args.rank and len(args.rank) > 1:
    args.rank = rank_code(rank)

seq_ids = defaultdict(list)
extract_ids = set()

# remember to use filtered nodes.dmp and names.dmp
def load_taxonomy(db_prefix):
    name_map = {}
    rank_map = {}
    child_lists = defaultdict(list)
    #read the taxonomy .dmp to and create or dict
    if not os.path.exists(db_prefix+"/taxonomy/name_map.json") or \
        not os.path.exists(db_prefix+"/taxonomy/rank_map.json") or \
        not os.path.exists(db_prefix+"/taxonomy/child_lists.json"):
        print ("Map files don't exist, creating json...", file=sys.stderr)
        with open(db_prefix+"/taxonomy/names_trimmed.dmp", 'r') as name_file:
            for line in name_file:
                node_id, name = line.strip().split('|')
                node_id = node_id.strip()
                name = name.strip()
                name_map[node_id] = name
        with open(db_prefix+"/taxonomy/nodes_trimmed.dmp", 'r') as nodes_file:
            for line in nodes_file:
                node_id, parent_id, rank = line.strip().split('|')
                node_id = node_id.strip()
                parent_id = parent_id.strip()
                rank = rank.strip()
                if node_id == '1':
                    parent_id = '0'
                child_lists[parent_id].append(node_id)
                rank_map[node_id] = rank
        #save our dicts as json
        with open(db_prefix+"/taxonomy/name_map.json",'w') as name_map_file, \
                open(db_prefix+"/taxonomy/rank_map.json",'w') as rank_map_file, \
                open(db_prefix+"/taxonomy/child_lists.json",'w') as child_lists_file:
            json.dump(name_map,name_map_file)
            json.dump(rank_map, rank_map_file)
            json.dump(child_lists,child_lists_file)
    else: #load the json
        with open(db_prefix+"/taxonomy/name_map.json",'r') as name_map_file, \
                open(db_prefix+"/taxonomy/rank_map.json",'r') as rank_map_file, \
                open(db_prefix+"/taxonomy/child_lists.json",'r') as child_lists_file:
            name_map = json.load(name_map_file)
            rank_map = json.load(rank_map_file)
            child_lists = json.load(child_lists_file)
        
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

def extract_fasta_from_id(fileout, id_list, seqfile):
    num_seq_to_extract = len(id_list)
    with open(fileout+".fa", 'w') as fout:
        for rec in SeqIO.parse(seqfile, 'fasta'):
            if rec.id in id_list: # as set is more efficient than a list
                #see https://wiki.python.org/moin/TimeComplexity
                num_seq_to_extract -= 1
                SeqIO.write(rec, fout,  'fasta')
                if num_seq_to_extract == 0:
                    break
        if num_seq_to_extract > 0:
            print ( "Warning, EOF reached but", num_seq_to_extract, "sequences remained", file=sys.stderr)

def extract_bam_from_id(fileout, id_list, seqfile):
    num_seq_to_extract = len(id_list)
    with AlignmentFile(seqfile, 'rb', check_sq=False) as bam_in, \
        AlignmentFile(fileout+".bam", 'wb', template=bam_in) as fout:
        for read in bam_in.fetch(until_eof=True):
            if read.query_name in id_list: # as set is more efficient than a list
                #see https://wiki.python.org/moin/TimeComplexity
                if not read.is_paired:
                    num_seq_to_extract -= 1
                elif read.is_read2: # decrease counter only if we see the second read of our pair
                    num_seq_to_extract -= 1
                fout.write(read)
                if not num_seq_to_extract:
                    break

def extract_seq_from_id(fileout, id_list, seqfile, data_type='bam'):
    if seqfile.endswith("fasta") or seqfile.endswith("fa") or seqfile.endswith("fas"):
        data_type = 'fasta'
    if data_type == 'fasta': extract_fasta_from_id(fileout, id_list, seqfile)
    elif data_type == 'bam': extract_bam_from_id(fileout, id_list, seqfile)

def dfs_report (node, depth):
    global extract_ids # we share this list through the recursive calls
    t_counts, c_counts, rank = taxo_counts[node], clade_counts[node], rank_map[node]
    if (not c_counts and not args.zeros):
        return
    c_counts_percent = round(c_counts * 100 / seq_count, 2)
    #filter on min seqences on clade
    #filter on min percent
    #filter on rank
    if (not args.rank or args.rank == rank_code(rank)) and (c_counts >= args.min and c_counts_percent >= args.minp):
        print ("{:6.2f}\t{}\t{}\t{}\t{}\t{}{}".format(
            c_counts_percent,
            c_counts,
            t_counts,
            rank_code(rank),
            node,
            "  " * depth,
            name_map[node]))
    # start saving the sequence mames for this clade
    if args.rank == rank_code(rank): extract_ids = set()
    children = child_lists.get(node,[])
    if len(children):
        sorted_children = sorted(children, key=lambda k: clade_counts[k], reverse=True)
        #format output only if not filtered by rank
        if not args.rank : depth += 1
        for child in sorted_children:
            dfs_report(child, depth)
    # we want to extract up to a certain clade from a certain rank,
    # if there is a min sequences to extract, and only if a ref file is provided
    if args.extractFile:
        if t_counts:# add only if the node has sequences assigned to it
            # a set is more efficient than a list: see https://wiki.python.org/moin/TimeComplexity
            extract_ids = extract_ids.union(seq_ids[node])
        if (node in args.clades or rank_code(rank) == args.rank) and \
            (c_counts_percent >= args.minp and len(extract_ids) >= args.min):
            print ("Extracting",len(extract_ids),"sequences for",name_map[node], file=sys.stderr)
            if "fa.kraken" in args.infile:
                suffix_length = len("fa.kraken")
            elif "fasta.kraken" in args.infile:
                suffix_length = len("fasta.kraken")
            else:
                suffix_length = len("kraken")
            # the names contains whitespaces
            extract_seq_from_id(args.infile[:-suffix_length]+name_map[node].replace(' ','_'), \
                                extract_ids, args.extractFile)
            extract_ids = set()

def dfs_summation(node):
    children = child_lists.get(node,[])
    if len(children):
        for child in children:
            dfs_summation(child)
            clade_counts[node] += clade_counts.get(child, 0)

name_map, rank_map, child_lists = load_taxonomy(db_prefix)

print("Map files loaded", file=sys.stderr)

seq_count = 0 # init the number of sequences
taxo_counts = defaultdict(int) # every new entry will be initialized to 0

with open(args.infile, 'r', newline='') as krakenfile:
    kfile = reader(krakenfile, delimiter='\t')
    for row in kfile:
        taxo_counts[row[2]] += 1
        seq_count += 1
        seq_ids[row[2]].append(row[1])

print(args.infile,"parsed", file=sys.stderr)

classified_count = seq_count - taxo_counts[0]
clade_counts = taxo_counts.copy()

if args.clades:
    for clade in args.clades:
        dfs_summation(clade)
else:
    dfs_summation('1')

unclassified_percent = 100
if seq_count:
    unclassified_percent = clade_counts.get(0) * 100 / seq_count
if  not args.clades and not args.rank:
    print ("{:6.2f}\t{}\t{}\t{}\t{}\t{}{}".format(
        unclassified_percent,
        clade_counts.get(0), taxo_counts[0],
        "U", 0, "", "unclassified"))

if args.clades:
    for clade in args.clades:
        dfs_report(clade, 0)
else:
    dfs_report('1', 0)
