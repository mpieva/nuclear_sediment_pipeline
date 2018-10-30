import gzip
from Bio import SeqIO, SeqRecord
import random
import sys
import pysam
import re
import argparse


# 30,000,000 with 99% bacteria => 99% bacteria, 0.5% hyena, 0.3% cow, 0.15% pig, 0.04% mammoth, 0.01% Human
# 0.5 -> 95912 0.3% -> 9913, 0.15% -> 9823, 0.04 -> 9785 0.01,  -> 9606

random.seed()


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
        sequence[-nbases+sequence[-nbases:].index("C")] = "T"
    return sequence.toseq()


# the idea is to generate a set of coordinates/size pairs and only take those substrings:
# sample 0:length_genome (N = #chunks desired)
# sample N sizes (N = chunks desired)
def chunk_fast(record, n_samples, vcf_in=None, chrom=None, individual=0, deaminate=0):
    try:
        positions = random.sample(range(0, len(record)-35), n_samples)
    except:
        print("sample to small", len(record)-35, n_samples, file=sys.stderr)
        positions = range(0, len(record)-35)
    # length = random.choices(range(35,100), k = n_samples) #only from python6 onwards
    length = [random.choice(range(35, 100)) for k in range(n_samples)]
    all_samples = []
    for pos, l in zip(positions, length):
        while 'N' in record[pos:pos+l]:
            # print ("Unknown base (N) in read at {} {}, resampling...".format(pos,l), file=sys.stderr) # debug purpose
            pos = random.choice(range(0, len(record)-35))
            l = random.choice(range(35, 100))
        # if not('N' in record[pos:pos+l]): # control for sequences without N's, deprecated by the 'while' statement above
        sample = record[pos:pos+l]
        if vcf_in:  # insert variation from VCF
            sequence = sample.seq.tomutable()
            for vcf_rec in vcf_in.fetch(chrom, start=pos, end=pos+l):
                # Insert alt only if GT != 0/0, use the first individual by default
                if any(vcf_rec.samples[0]['GT']):
                    if len(vcf_rec.alts[0]) == 1:  # only SNPs
                        sequence[vcf_rec.start-pos] = vcf_rec.alts[0]
                # print(record[vcf_rec.start-3:vcf_rec.start+4].seq,vcf_rec.ref,vcf_rec.alts,vcf_rec.start, pos, l) # debug purpose
            sample.seq = sequence.toseq()
        if deaminate:
            sample.seq = simulate_deamination(
                sample.seq.tomutable(), deaminate)
        all_samples += [(sample, pos)]
    return all_samples


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
                        int(s/full_size * num_seq)+1 for s in fa.lengths[:n_chromosomes]]
                else:
                    full_size = sum(fa.lengths)
                    size_percent = [
                        int(s/full_size * num_seq)+1 for s in fa.lengths]
            return (size_percent)
    except OSError:  # fasta is not bgzip'd do a naive fallback
        print("Error, fasta file not indexed, doing naive sampling...", file=sys.stderr)
        if n_chromosomes:
            return [int(num_seq / n_chromosomes)+1]*n_chromosomes
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
    parser.add_argument(
        '--vcf',  help='VCF containing the variation to be included in the samples')
    parser.add_argument('file_in', metavar="genome.fa.gz, bgzip'd and indexed")
    parser.add_argument('--chromosomes', type=int,
                        help='How many chromosomes are expected', default=0)
    parser.add_argument('--deaminate', type=int,
                        help='How many bases should be deaminated of each end', default=0)
    parser.add_argument('--individual', type=int,
                        help='The individual used for variation', default=0)
    args = parser.parse_args()

    print("Loading dataset...", end="", file=sys.stderr)

    num_reads_to_sample = estimate_read_distribution(
        args.file_in, args.num_seq, args.chromosomes)
    with gzip.open(args.file_in, "rt") as file_in:
        all_chunks = []
        # if args.chromosomes: # we iterate record by record
        #    list_records = SeqIO.parse(file_in, "fasta")
        # over sampling here because some samples will be N only and thus discarded
        # TODO, to naive, chromosome should be sampled with respect to their size
        #    num_reads_to_sample = int(args.num_seq  / (max(1, args.chromosomes-1)))
        #
        #    print("Done.\nWill sample an average of {} reads per sequence.".format(num_reads_to_sample), file=sys.stderr)
        # else:
        #    list_records = list(SeqIO.parse(file_in, "fasta"))
        # TODO: sample % relative to the record size / total size (more reads from chr1 than chr21 for example
        #    num_reads_to_sample = int(args.num_seq  / (len(list_records)-1))
        #    print("Done\n{} sequences loaded.\nWill sample an average of {} reads per sequence".format(len(list_records), num_reads_to_sample), file=sys.stderr)
        # for n_record, record in enumerate(list_records, 1):
        # TODO, check if specifying the chunk size can be directly in map...nope, size is always 1
        vcf_in = chrom = None
        if args.vcf:
            vcf_in = pysam.VariantFile(args.vcf)
            # subset the VCF to 1 individual only
            vcf_in.subset_samples([vcf_in.header.samples[args.individual]])
            # our chromosome can contain up to 3 characters (most likely 2 max: e.g. chromosome 21,)
            p = re.compile("chromosome \w{1,3},", re.IGNORECASE)
        # list_records):
        for num_record, record in enumerate(SeqIO.parse(file_in, "fasta")):
            # we want reads only to assigned chromosomes
            if args.chromosomes:
                print("Parsing chromosome", num_record,
                      end="\r", file=sys.stderr)
                if num_record > args.chromosomes:
                    # assume the genome is sorted, with chromosomes first...
                    break
            if args.vcf:
                chrom = p.search(record.description)
                if not chrom:
                    chrom = record.id  # in case of scaffold genome
                else:
                    chrom = chrom.group()[len('chromosome '):-1]
            all_chunks += chunk_fast(
                record, num_reads_to_sample[num_record], vcf_in, chrom, deaminate=args.deaminate)
            if num_record % 100 == 99:  # show progress
                print(num_record+1, "sequences parsed...",
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
