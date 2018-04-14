from Bio import SeqIO
import random
import sys
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Split a genome into chuncks, add mutations')
parser.add_argument('--mut_rate', default=0.001, type=float,
                    help='Mutation rate')
parser.add_argument('--num_seq', default=0, type=int,
                    help='Number of sequences to output')
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


def main():
    all_chunks = []
    #first_record = next(SeqIO.parse("/mnt/sequencedb/Refseq/viral_genomes.fa", "fasta"))
    print("Searching for specie", args.specie, file=sys.stderr) 
    for record in SeqIO.parse(args.file_in, "fasta"):
        if check_kraken_id(record, args.specie): #start thread
            print("Found sequence for specie", args.specie, record.id, file=sys.stderr) 
            i=0
            tmp_chunks = []
            while i < len(record):
                # discard sequences with N only
                for nuc in record[i:]:
                    if nuc == 'N': i+=1
                    else: break
                #default size is 35<>100
                chunk_size = get_random_length()
                #print ("Chunk",i,i+chunk_size, record[i:i+chunk_size].seq)
                #get a subseq of the record
                #chunk_record = record[i:i+chunk_size]                
                #all_chunks.append([chunk_record, i])
                
                #reduce the size of the array for memory efficiency
                tmp_chunks.append([record[i:i+chunk_size], i])
                i += chunk_size
            if args.num_seq:
                all_chunks += random.sample(tmp_chunks, args.num_seq)
            else:
                all_chunks += tmp_chunks
            break
        # assume sequence are consecutive for each species 
        elif len(all_chunks):#join threads
            break
    if args.num_seq:
        all_chunks = random.sample(all_chunks, args.num_seq)
     
    with open(args.outfile, 'w') as file_out:
        for chunk_record, pos in all_chunks:
            # add mutations to sequence only now, as we work on a subset (runs faster?)
            if args.mut_rate:
                chunk_record.seq = insert_mutations(chunk_record.seq)
            chunk_record.id = "{}|{}_{}".format(chunk_record.id, pos, len(chunk_record))
            SeqIO.write(chunk_record, file_out, "fasta")
            break

if __name__ == "__main__":
    main()
