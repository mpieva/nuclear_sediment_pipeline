"""Parse a kraken output file and generate a report and possibly extract reads for selected clades. (Adapted from original kraken-report.pl)
"""

import sys
import gzip
from csv import reader
from Bio import SeqIO
from pysam import AlignmentFile
from collections import defaultdict
import argparse
from pathlib import Path
import os
import json
import contextlib
#grep 'scientific name' names.dmp |cut -d'|' -f 1,2 |gzip -c >names_trimmed.dmp
#cut -d '|' -f 1,2,3 nodes.dmp|gzip -c >nodes_trimmed.dmp

# remember to use filtered nodes.dmp and names.dmp
def load_taxonomy(db_prefix):
    """Create/Read a taxonomy maps into dicts
    """
    global name_map
    name_map = {}
    global rank_map
    rank_map = {}
    global child_lists
    child_lists = defaultdict(list)
    global name_clade_map
    parent_map = {}
    #read the taxonomy .dmp to and create or dict
    if not os.path.exists(db_prefix+"/taxonomy/name_map.json") or \
        not os.path.exists(db_prefix+"/taxonomy/rank_map.json") or \
        not os.path.exists(db_prefix+"/taxonomy/child_lists.json") or \
        not os.path.exists(db_prefix+"/taxonomy/parent_map.json"):
        print ("Map files don't exist, creating json...", file=sys.stderr)
        with gzip.open(db_prefix+"/taxonomy/names_trimmed.dmp.gz", 'rt') as name_file:
            for line in name_file:
                node_id, name = line.strip().split('|')
                node_id = node_id.strip()
                name = name.strip()
                name_map[node_id] = name
        with gzip.open(db_prefix+"/taxonomy/nodes_trimmed.dmp.gz", 'rt') as nodes_file:
            for line in nodes_file:
                node_id, parent_id, rank = line.strip().split('|')
                node_id = node_id.strip()
                parent_id = parent_id.strip()
                rank = rank.strip()
                if node_id == '1':
                    parent_id = '0'
                child_lists[parent_id].append(node_id)
                rank_map[node_id] = rank
                parent_map[node_id] = parent_id
        #save our dicts as json
        with open(db_prefix+"/taxonomy/name_map.json",'w') as name_map_file, \
                open(db_prefix+"/taxonomy/rank_map.json",'w') as rank_map_file, \
                open(db_prefix+"/taxonomy/child_lists.json",'w') as child_lists_file, \
                open(db_prefix+"/taxonomy/parent_map.json",'w') as parent_map_file:
            json.dump(name_map,name_map_file)
            json.dump(rank_map, rank_map_file)
            json.dump(child_lists,child_lists_file)
            json.dump(parent_map, parent_map_file)
    else: #load the json
        with open(db_prefix+"/taxonomy/name_map.json",'r') as name_map_file, \
                open(db_prefix+"/taxonomy/rank_map.json",'r') as rank_map_file, \
                open(db_prefix+"/taxonomy/child_lists.json",'r') as child_lists_file:
            name_map = json.load(name_map_file)
            rank_map = json.load(rank_map_file)
            child_lists = json.load(child_lists_file)
    name_clade_map = {v: k for k, v in name_map.items()}
    #return (name_map, rank_map, child_lists, name_clade_map)

def rank_code(rank):
    """Translate ranks into single letters code
    """
    if rank == "species": return "S"
    if rank == "genus": return "G"
    if rank == "family": return "F"
    if rank == "order": return "O"
    if rank == "class": return "C"
    if rank == "phylum": return "P"
    if rank == "kingdom": return "K"
    if rank == "superkingdom": return "D"
    return "-"

def get_taxonomy_str(taxid):
    """Generate the full taxonomy from a specific clade
    
    Parameters
    ----------
    taxid: str
    
    Returns
    -------
    str
    """
    taxid_string = known_taxonomy_strings.get(taxid, False)
    if not taxid_string:
        nodes = []
        while (taxid != '0'):
            nodes += [name_map[taxid]]
            taxid = parent_map[taxid]
        taxid_string = ';'.join(nodes[::-1])
        known_taxonomy_strings[taxid] = taxid_string
    return taxid_string
    

@contextlib.contextmanager
def _ret_file(f):
    yield f
    
def extract_fasta_from_id(fileout, id_list, seqfile, min_length):
    """Extract reads assigned to specific taxa.
    
    Parameters
    ----------
    fileout: str
        Filename to write into
    id_list: list of 
    """
    if seqfile.endswith('a') or seqfile.endswith('a.gz'):
        file_type = "fasta"
        file_suffix = '.fa'
    elif seqfile.endswith('q') or seqfile.endswith('q.gz'):
        file_type = "fastq"
        file_suffix = '.fq'
    with open(fileout+file_suffix, 'w') as fout, \
    gzip.open(seqfile, "rt") if seqfile.endswith("gz") else _ret_file(seqfile) as seqfile:
        # working with a generator expression, may be better memory-wise
        input_seq_iterator = SeqIO.parse(seqfile, file_type)
        fasta_seq_iterator = (rec for rec in input_seq_iterator if rec.id in id_list and len(rec) >= min_length)
        count = SeqIO.write(fasta_seq_iterator, fout, file_type)
    if len(id_list) != count: # sanity check you may want to extract from a demultiplexed file
        print("Warning, EOF reached but", len(id_list) - count, "sequences remained, is extractFile the original source?", file=sys.stderr)

def extract_bam_from_id(fileout, id_list, seqfile, min_length):
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
                if read.query_length >= min_length:
                    fout.write(read)
                if not num_seq_to_extract:
                    break

def extract_seq_from_id(fileout, id_list, seqfile, data_type='bam', min_length=0):
    if seqfile.endswith("fasta") or seqfile.endswith("fa") or seqfile.endswith("fas") or seqfile.endswith("fq") or seqfile.endswith("fastq") or seqfile.endswith("gz"):
        data_type = 'fasta'
    if data_type == 'fasta': extract_fasta_from_id(fileout, id_list, seqfile, min_length)
    elif data_type == 'bam': extract_bam_from_id(fileout, id_list, seqfile, min_length)

def dfs_report (node, depth, related=[], infile="", extractFile=None, outdir="", zeros=False, clades=[], target_rank=None, min_count=0, minp=0.0, min_length=0):
    global extract_ids # we share this list through the recursive calls
    t_counts, c_counts, rank = taxo_counts[node], clade_counts[node], rank_map[node]
    if (not c_counts and not zeros):
        return
    c_counts_percent = round(c_counts * 100 / seq_count, 2)
    #filter on min seqences on clade
    #filter on min percent
    #filter on rank
    if (not target_rank or target_rank == rank_code(rank)) and (c_counts >= min_count and c_counts_percent >= minp):
        if node not in related: # TODO not in excluded, implement an 'excluded' switch
            print ("{:6.2f}\t{}\t{}\t{}\t{}\t{}{}".format(
                c_counts_percent,
                c_counts,
                t_counts,
                rank_code(rank),
                node,
                "  " * depth,
                name_map[node]))
    # start saving the sequence mames for this clade
    if target_rank == rank_code(rank): extract_ids = set()
    children = child_lists.get(node,[])
    if len(children):
        sorted_children = sorted(children, key=lambda k: clade_counts[k], reverse=True)
        #format output only if not filtered by rank
        if not target_rank : depth += 1
        for child in sorted_children:
            #dfs_report(child, depth)
            dfs_report(child, depth, related, infile, extractFile, outdir, zeros, clades, target_rank, min_count, minp, min_length)
    # we want to extract up to a certain clade from a certain rank,
    # if there is a min sequences to extract, and only if a ref file is provided
    if extractFile:   
        outdir = Path(outdir) 
        if not outdir.exists():
            outdir.mkdir(parents=True)
        if t_counts:# add only if the node has sequences assigned to it
            # a set is more efficient than a list: see https://wiki.python.org/moin/TimeComplexity
            extract_ids = extract_ids.union(seq_ids[node])
        if (node in clades or rank_code(rank) == target_rank) and \
            (c_counts_percent >= minp and len(extract_ids) >= min_count):
            print ("Extracting",len(extract_ids),"sequences for",name_map[node], file=sys.stderr)
            if "fa.kraken" in infile or "fq.kraken" in infile:
                suffix_length = len("fa.kraken")
            elif "fasta.kraken" in infile or "fastq.kraken" in infile:
                suffix_length = len("fasta.kraken")
            elif "bam.kraken" in infile:
                suffix_length = len("bam.kraken")
            else:
                suffix_length = len("kraken")
            # the names contains whitespaces
            extract_seq_from_id(str(outdir / Path(infile).name[:-suffix_length])+name_map[node].replace(' ','_'), \
                                extract_ids, extractFile, min_length)
            extract_ids = set()

def dfs_summation(node):
    children = child_lists.get(node,[])
    if len(children):
        for child in children:
            dfs_summation(child)
            clade_counts[node] += clade_counts.get(child, 0)

#this function will discard child clades in order to have a proper summation
def dfs_related(node, node_list):
    res = []
    children = child_lists.get(node,[])
    if len(children):
        #iterate through all children
        for child in children:
            if child in node_list: res+=[child]
            # recursively look for children
            res += dfs_related(child, node_list)
    return res
    
def _main():
    parser = argparse.ArgumentParser(description='Create a report from a kraken output. Optionally extract reads')
    parser.add_argument('--db', required=True,
                        help='The kraken database to use')
    parser.add_argument('--zeros', action='store_true',
                        help='Show also 0')
    parser.add_argument('--clades', default=[],
                        help='Select only specified clades (comma separated)')
    parser.add_argument('--minp', default=0.0, type=float,
                        help='Filter on the minimum percent of sequences for this clade')
    parser.add_argument('--min_count', default=0, type=int,
                        help='Filter on the minimum sequences for this clade')
    parser.add_argument('--rank',  help='Only return clades for specified rank')
    parser.add_argument('--translate',  help='Output for "translate" (read -> lineage)')
    parser.add_argument('--extractFile',  help='File where to extract sequence from')
    parser.add_argument('--min_length', default=0, type=int, help='Minimum length filter')
    parser.add_argument('infile', metavar="kraken.output")
    parser.add_argument('--outdir', default="", help='Extracted reads directory')

    args = parser.parse_args()

    db_prefix = os.path.abspath(args.db)
    if args.rank and len(args.rank) > 1:
        args.rank = rank_code(rank)
    
    global seq_ids
    seq_ids = defaultdict(list)
    #extract_ids = set()
    load_taxonomy(db_prefix)
    #name_map, rank_map, child_lists, node_name_map = load_taxonomy(db_prefix)
    known_taxonomy_strings = {}
    if args.translate:
        with open(db_prefix+"/taxonomy/parent_map.json",'r') as parent_map_file:
            parent_map = json.load(parent_map_file)

    print("Map files loaded", file=sys.stderr)

    if args.clades: # handle providing multiple clades, comma separated
        args.clades = args.clades.split(",")
        for idx, clade in enumerate(args.clades): # translate taxa names to number
            if not clade in name_map:
                try:            
                    args.clades[idx] = node_name_map[clade.replace("_"," ")]
                except KeyError:
                    print("Specified taxa {} not found, exiting", file=sys.stderr)
                    exit(1)

        args.clades = set(args.clades)
    
    global seq_count
    seq_count = 0 # init the number of sequences
    global taxo_counts
    taxo_counts = defaultdict(int) # every new entry will be initialized to 0

    with open(args.infile, 'r', newline='') as krakenfile, \
        open(args.translate, "w") if args.translate else _ret_file(None) as translate:
        kfile = reader(krakenfile, delimiter='\t')
        for row in kfile:
            taxo_counts[row[2]] += 1
            seq_count += 1
            seq_ids[row[2]].append(row[1])
            if args.translate and row[0].startswith('C'):
                    print (row[1], get_taxonomy_str(row[2]), sep="\t", file=translate)

    print(args.infile,"parsed", file=sys.stderr)

    classified_count = seq_count - taxo_counts[0]
    global clade_counts
    clade_counts = taxo_counts.copy()

    if args.clades:
        #do the summation only once for each clade,
        # that means, if we specify clades related to each other:
        # e.g. 9606 9605, only the higher clade will be used
        # as the descendant one will be recursively computed
        related_clades=set()
        for node in args.clades:
            related_clades = related_clades.union(dfs_related(node, args.clades))
        unrelated_clades = args.clades.difference(related_clades)
        for clade in unrelated_clades:
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
        related_clades=set()
        for node in args.clades:
            related_clades = related_clades.union(dfs_related(node, args.clades))
        for clade in args.clades:
        #dfs_report (node, depth, related=[], infile="", extractFile=None, outdir="", zeros=False, clades=[], target_rank=None, min_count=0, minp=0.0, min_length=0)
            dfs_report(clade, 0, related_clades, args.infile, args.extractFile, args.outdir, args.zeros, args.clades, args.rank, args.min_count, args.minp, args.min_length)
    else:
        dfs_report('1', 0, [], args.infile, args.extractFile, args.outdir, args.zeros, args.clades, args.rank, args.min_count, args.minp, args.min_length)

if __name__ == "__main__":
    name_map = rank_map = child_lists = node_name_map = clade_counts = taxo_counts = seq_count = extract_ids = seq_ids = None
    _main()
