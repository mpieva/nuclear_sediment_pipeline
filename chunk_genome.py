from Bio import SeqIO, SeqRecord
import random
import sys
import numpy as np
from multiprocessing import Pool
import argparse

parser = argparse.ArgumentParser(description='Split a genome into chuncks, add mutations')
parser.add_argument('--mut_rate', default=1/10000, type=float,
                    help='Mutation rate')
parser.add_argument('--num_seq', default=0, type=int,
                    help='Number of sequences to output')
parser.add_argument('--nthreads', default=4, type=int,
                    help='Number of threads to use')
parser.add_argument('--outfile',  help='Fasta file where to extract sequence')
parser.add_argument('--specie',  help='Specie Taxa')
parser.add_argument('file_in', metavar="genome.fa")

args = parser.parse_args()


# 3000,000,000 with 99% bacteria => 99% bacteria, 0.5% hyena, 0.3% cow, 0.15% pig, 0.04% mammoth, 0.01% Human
# 0.5 -> 95912 0.3% -> 9913, 0.15% -> 9823, 0.04 -> 9785 0.01,  -> 9606

alphabet = set(["A","C","G","T"])
random.seed()

def check_kraken_id (seq, kraken_id='9606'):
    return seq.id.split("|")[2] == kraken_id

def get_random_length(min_length = 35, max_length = 100):
    return random.randint(min_length, max_length)
    
#todo insert C->T method
def insert_mutations(seq):
    s = seq.tomutable()
    # creates an array of position with uniform [0,1) values
    pos_mut_threshold = np.random.rand(len(s))
    # if the element in the array is above threshold, insert a mutation in the sequence
    for pos in np.where(pos_mut_threshold <= args.mut_rate)[0]:
        # we need pos.item() otherwise s[pos] returns a new sequence
        s[pos.item()] = random.sample(alphabet - set(s[pos]), 1)[0]
    return s

def iterate_through_record(record):
    i=0
    tmp_chunks = []
    while i < len(record):
        # discard sequences with N only
        if record[i] == 'N':
            i+=1
            continue
        #default size is 35<>100
        chunk_size = get_random_length()
        
        #reduce the size of the array for memory efficiency
        tmp_chunks.append([record[i:i+chunk_size], i])
        i += chunk_size
    
    # control if our chunk contains enough reads
    try:
        if args.num_seq:
            return random.sample(tmp_chunks, int(args.num_seq/args.nthreads)+1)
    except ValueError:
        pass
    return tmp_chunks

def main():
    all_chunks = []
    #first_record = next(SeqIO.parse("/mnt/sequencedb/Refseq/viral_genomes.fa", "fasta"))
    print("Searching for specie", args.specie, file=sys.stderr) 
    for record in SeqIO.parse(args.file_in, "fasta"):
        if check_kraken_id(record, args.specie): #start thread
            print("Found sequence for taxa", args.specie, record.id, file=sys.stderr) 
            with Pool(args.nthreads) as p:
                #TODO, check if specifying the chunk size can be directly in map...nope, size is always 1
                thread_chunk_size = int(len(record) / args.nthreads)+1
                tmp_chunks = p.map(iterate_through_record, [record[i*thread_chunk_size:(i+1)*thread_chunk_size] for i in range(args.nthreads)])
                #TODO sample x chunks directly (random [start positions and length times] * num_seq)
            for thread_chunks in tmp_chunks:
                all_chunks += thread_chunks
        # assume sequence are consecutive for each species 
        elif len(all_chunks):
            break
    print("All sequences parsed, merging...", file=sys.stderr)
    try:
        if args.num_seq:
            all_chunks = random.sample(all_chunks, args.num_seq)
    except ValueError:
        print("Returning only {} sequences".format(len(all_chunks)),file=sys.stderr)
    with open(args.outfile, 'w') as file_out:
        if args.mut_rate:
            record_it = (SeqRecord.SeqRecord(insert_mutations(record.seq), id="{}|{}_{}".format(record.id, pos, len(record)), description = record.description) for record, pos in all_chunks)
        else:
            record_it = (SeqRecord.SeqRecord(record.seq, id="{}|{}_{}".format(record.id, pos, len(record)), description = record.description) for record, pos in all_chunks)
        SeqIO.write(record_it, file_out, "fasta")
if __name__ == "__main__":
    main()
