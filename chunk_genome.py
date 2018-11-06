import gzip
from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import random
import sys
import pysam
import re
import argparse


# 30,000,000 with 99% bacteria => 99% bacteria, 0.5% hyena, 0.3% cow, 0.15% pig, 0.04% mammoth, 0.01% Human
# 0.5 -> 95912 0.3% -> 9913, 0.15% -> 9823, 0.04 -> 9785 0.01,  -> 9606

# tabl['A'].loc[30][tabl['A'].loc[30]>=0.0003].idxmin()

random.seed()


def get_sequence_with_substitution(sequence):
    read_length = len(sequence)
    choices = np.random.random(read_length)
    positions = np.concatenate([np.arange(15), np.arange(29, 14, -1)])
    positions = np.insert(positions, 15, np.full(read_length - 30, 30))
    newSeq = [''] * read_length  # working on list is faster
    #same_nuc = 1-min([min([tabl[x][pos][x] for pos in range(1,30)]) for x in tabl])
    #same_nuc = {nc: 1-prob for nc, prob in zip(['A','C','G','T'],min([[tabl[x][pos][x] for x in tabl] for pos in range(1,30)]))}

    # max_other = min([tabl[x][tabl[x]<same_nuc].max() for x in tabl])# this
    # is the highest value. Anything above we sample the same nuc
    for idx, nc in enumerate(sequence):
        pos = positions[idx]
        choice = choices[idx]
        # TODO optimize speed for pos '30': if prob > max(30), no matter which
        # base, just use nc
        if choice > same_nuc[nc]:
            newSeq[idx] = nc
        else:
            try:
                newSeq[idx] = tabl[nc].loc[pos][
                    tabl[nc].loc[pos] <= choice].idxmax()
            except:
                newSeq[idx] = tabl[nc].loc[pos].idxmin()
                #newSeq[idx] = nc
    return Seq(''.join(newSeq))
    # returns the base relative to our random number for position pos given
    # the real base is nc

# specific code to work with mappability track...


def get_map_pos(n_samples, map_file="/tmp/fred/map_track.bed.gz"):
    with pysam.TabixFile(map_file, parser=pysam.asBed()) as tbx:
        reads = random.sample([row for row in tbx.fetch() if (
            int(row[2]) - int(row[1]) < 100) and (int(row[2]) - int(row[1]) > 35)], n_samples)
    list_records = SeqIO.parse(
        '/mnt/solexa/Genomes/hg19_1000g/whole_genome.fa', "fasta")
    chromosomes = 24  # 22 + X + Y
    all_chunks = []
    for num_record, record in enumerate(list_records, 1):
        if num_record > chromosomes:
            break
        chrom = record.id
        chrom_sample = [row for row in reads if row[0] == chrom]
        print(chrom, len(chrom_sample), file=sys.stderr)
        for pos in chrom_sample:
            sample = record[int(pos[1]):int(pos[2])]
            all_chunks += [(sample, pos[1])]

    record_it = (SeqRecord.SeqRecord(record.seq, id="{}|{}_{}|SGDP".format(record.id, pos, len(record)),
                                     description=" ".join(record.description.split(' ')[1:])) for record, pos in all_chunks)
    with open('/tmp/fred/human.fa', 'w') as file_out:
        SeqIO.write(record_it, file_out, "fasta")


def simulate_deamination(sequence, nbases=3):
    while "C" in sequence[:nbases]:  # 5' C to T
        sequence[sequence.index("C")] = "T"
    while "C" in sequence[-nbases:]:  # 3' C to T
        sequence[-nbases + sequence[-nbases:].index("C")] = "T"
    return sequence.toseq()


def mutate_unif(sequence, unif):
    same_nuc = set("ACGT")
    read_length = len(sequence)
    # choices is an array of True/False
    # we will mutate if or sample is < unif
    choices = np.random.random(read_length) < unif
    newSeq = [''] * read_length
    for idx, nc in enumerate(sequence):
        mutate = choices[idx]
        if mutate:
            # we call tuple in order to do random.choice (and remove the current nc from the possibilities)
            newSeq[idx] = random.choice(tuple(same_nuc.difference(nc)))
        else:
            newSeq[idx] = nc
    return Seq(''.join(newSeq))


# the idea is to generate a set of coordinates/size pairs and only take those substrings:
# sample 0:length_genome (N = #chunks desired)
# sample N sizes (N = chunks desired)
def chunk_fast(record, n_samples, vcf_in=None, chrom=None, individual=0, unif=False, len_distrib=False, deaminate=0, minlength=35, maxlength=100):
    try:
        positions = random.sample(range(0, len(record)-minlength), n_samples)
    except:
        print("sample too small {}bp, sampling {} reads with replacements".format(len(record) -
                                                                                  minlength, n_samples), file=sys.stderr)
        positions = [random.choice(range(0, len(record) - minlength))
                     for _ in range(n_samples)]
    if len_distrib:  # we gave a file with distribution per length
        p = read_len_distrib(len_distrib, minlength, maxlength)
        length = np.random.choice(
            np.arange(minlength, maxlength+1), n_samples, p=p)
    else:
        # length = random.choices(range(minlength, maxlength), k = n_samples) #only from python6
        # onwards
        length = [random.choice(range(minlength, maxlength+1))
                  for k in range(n_samples)]
    all_samples = []
    for pos, l in zip(positions, length):
        while 'N' in record[pos:pos + l]:
            # print ("Unknown base (N) in read at {} {},
            # resampling...".format(pos,l), file=sys.stderr) # debug purpose
            pos = random.choice(range(0, len(record) - minlength))
            if len_distrib:
                l = np.random.choice(
                    np.arange(minlength, maxlength+1), 1, p=p)[0]
            else:
                l = random.choice(range(minlength, maxlength+1))
        # if not('N' in record[pos:pos+l]): # control for sequences without
        # N's, deprecated by the 'while' statement above
        sample = record[pos:pos + l]
        if vcf_in:  # insert variation from VCF
            sequence = sample.seq.tomutable()
            if isinstance(vcf_in, pysam.VariantFile):  # we use a VCF
                for vcf_rec in vcf_in.fetch(chrom, start=pos, end=pos + l):
                    # Insert alt only if GT != 0/0, use the first individual by
                    # default
                    if any(vcf_rec.samples[0]['GT']):
                        if len(vcf_rec.alts[0]) == 1:  # only SNPs
                            sequence[vcf_rec.start - pos] = vcf_rec.alts[0]
                sample.seq = sequence.toseq()
                # print(record[vcf_rec.start-3:vcf_rec.start+4].seq,vcf_rec.ref,vcf_rec.alts,vcf_rec.start,
                # pos, l) # debug purpose
            else:  # we use a variation probability matrix
                read_length = len(sequence)
                sample.seq = get_sequence_with_substitution(sequence)
        if deaminate:
            sample.seq = simulate_deamination(
                sample.seq.tomutable(), deaminate)
        if unif:
            sample.seq = mutate_unif(sample.seq.tomutable(), unif)
        all_samples += [(sample, pos)]
    return all_samples


def read_len_distrib(filename, minlen=35, maxlen=100):
    from csv import reader
    with open(filename, 'r', newline='') as csvfile:
        csvreader = reader(csvfile, delimiter='\t')
        first_row = csvreader.__next__()
        read_distrib = [0] * int(first_row[1])+[int(first_row[0])]
        read_distrib.extend([int(row[0]) for row in csvreader])
    if len(read_distrib) < maxlen:
        print("Read length distribution provides values only until",
              len(read_distrib), "; expected:", maxlen)
    s = sum(read_distrib[minlen:maxlen+1])
    return [r / s for r in read_distrib[minlen:maxlen+1]]


def estimate_read_distribution(file_in, num_seq, n_chromosomes=None):
    try:
        with pysam.FastaFile(file_in) as fa:
            print("The file contains", fa.nreferences,
                  "sequences", file=sys.stderr)
            if num_seq:
                if n_chromosomes:
                    print("Working only with the first",
                          n_chromosomes, file=sys.stderr)
                    full_size = sum(fa.lengths[:n_chromosomes])
                    size_percent = [
                        int(s / full_size * num_seq) + 1 for s in fa.lengths[:n_chromosomes]]
                else:
                    full_size = sum(fa.lengths)
                    size_percent = [
                        int(s / full_size * num_seq) + 1 for s in fa.lengths]
            return (size_percent)
    except OSError:  # fasta is not bgzip'd do a naive fallback
        print("Error, fasta file not indexed, doing naive sampling...",
              file=sys.stderr)
        if n_chromosomes:
            return [int(num_seq / n_chromosomes) + 1] * n_chromosomes
        else:
            print(
                "Naive sampling needs the option --chromosomes...Aborting...", file=sys.stderr)
            sys.exit(1)


def main():
    # Process command line
    parser = argparse.ArgumentParser(description='Split a genome into chuncks')
    parser.add_argument('--num_seq', default=0, type=int,
                        help='Number of sequences to output')
    parser.add_argument(
        '--outfile',  help='Fasta file where to extract sequence')
    parser.add_argument('--specie',  help='Specie Taxa')
    parser.add_argument('file_in', metavar="genome.fa.gz, bgzip'd and indexed")
    parser.add_argument('--chromosomes', type=int,
                        help='How many chromosomes are expected', default=0)
    parser.add_argument('--minlen', type=int,
                        help='Minimum read length (lengths are sampled uniformly)', default=35)
    parser.add_argument('--maxlen', type=int,
                        help='Maximum read length (lengths are sampled uniformly)', default=100)
    # TODO, control for option dependency properly
    parser.add_argument('--deaminate', type=int,
                        help='How many bases should be deaminated of each end', default=0)
    parser.add_argument('--unif', type=float,
                        help='Add mutations uniformly distributed', default=0.0)
    parser.add_argument('--individual', type=int,
                        help='The individual used for variation in you VCF', default=0)
    parser.add_argument(
        '--length', help='2 columns (#reads, length) TSV file containing read length distribution', default=0.0)
    # either provide a VCF OR a substitution matrix, you probably don't want
    # to deaminate if you choose the latter
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        '--vcf',  help='VCF containing the variation to be included in the samples')
    group.add_argument("--substitution_file", type=open,
                       help='File containing the substitution probability matrix')
    args = parser.parse_args()

    print("Loading dataset...", end="", file=sys.stderr)

    if args.substitution_file:
        global tabl
        t = pd.read_table(args.substitution_file, sep="\t",
                          header=None, names=["profile", "from", "to", "prob"])
        tabl = pd.pivot_table(t, values='prob', index=[
                              'profile', 'to'], columns=['from'], aggfunc=np.sum)
        # contains, for each nc, the highest probability for which we can have
        # a mutation
        global same_nuc
        # e.g.
        # 'A': 0.00797800000000004, -> ~0.79% chances to have a mutation, any sampled number above this would result in keeping the nc
        # this allows us to optimize the algorithm by not using the lookup
        # table and directly keep the nc.
        same_nuc = {nc: 1 - prob for nc, prob in zip(tabl.columns, min(
            [[tabl[x][pos][x] for x in tabl] for pos in range(1, 30)]))}
    num_reads_to_sample = estimate_read_distribution(
        args.file_in, args.num_seq, args.chromosomes)
    with gzip.open(args.file_in, "rt") as file_in:
        all_chunks = []
        vcf_in = args.substitution_file  # will be none if we use a VCF file
        chrom = None
        if args.vcf:
            vcf_in = pysam.VariantFile(args.vcf)
            # subset the VCF to 1 individual only
            vcf_in.subset_samples([vcf_in.header.samples[args.individual]])
            # our chromosome can contain up to 3 characters (most likely 2 max:
            # e.g. chromosome 21,)
        # used only when adding variation
        p = re.compile("chromosome \w{1,3},", re.IGNORECASE)
        # list_records):
        for num_record, record in enumerate(SeqIO.parse(file_in, "fasta")):
            # we want reads only to assigned chromosomes
            if args.chromosomes:
                print("Parsing chromosome", num_record +
                      1, end="\r", file=sys.stderr)
                if num_record > args.chromosomes - 1:  # num_record is 0-based
                    # assume the genome is sorted, with chromosomes first...
                    break
            if args.vcf or args.substitution_file:
                chrom = p.search(record.description)
                if not chrom:
                    chrom = record.id  # in case of scaffold genome
                else:
                    chrom = chrom.group()[len('chromosome '):-1]
                all_chunks += chunk_fast(record, num_reads_to_sample[
                                         num_record], vcf_in, chrom, unif=args.unif, deaminate=args.deaminate, len_distrib=args.length, minlength=args.minlen, maxlength=args.maxlen)
            if num_record % 100 == 99:  # show progress
                print(num_record + 1, "sequences parsed...",
                      end="\r", file=sys.stderr)
        if args.vcf:
            vcf_in.close()
        print("Sampling done, merging...", end="", file=sys.stderr)
        try:
            if len(all_chunks) > args.num_seq:
                all_chunks = random.sample(all_chunks, args.num_seq)
        except ValueError:
            print("Returning only {} sequences for {}".
                  format(len(all_chunks, record.id)), file=sys.stderr)
    print("Done\nWritting down records...", end="", file=sys.stderr)
    # create a new header which includes the read pos, read length
    record_it = (SeqRecord.SeqRecord(record.seq, id="{}|{}_{}".format(record.id, pos, len(record)),
                                     description=" ".join(record.description.split(' ')[1:])) for record, pos in all_chunks)
    if args.outfile:
        with open(args.outfile, 'w') as file_out:
            SeqIO.write(record_it, file_out, "fasta")
    else:
        SeqIO.write(record_it, sys.stdout, "fasta")
    print("Done", file=sys.stderr)


if __name__ == "__main__":
    main()
