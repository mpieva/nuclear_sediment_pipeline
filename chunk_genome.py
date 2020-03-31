"""Set of functions used to generate reads derived from a set of sequences (i.e. full genome FASTA file).

* You can specify a mutation rate, a VCF file in order to add known SNPs, or a mutation matrix obtained from SNPad.
* You can specify a length distribution file otherwise reads length will be uniformly sampled.
* Works is processed in parallel in order to decrease computation time.

"""

import gzip
from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq
import multiprocessing
# for creating partial function for feeding pool.map
from functools import partial
import numpy as np
import pandas as pd
import random
import sys
import pysam
import re
import argparse
import contextlib

# 30,000,000 with 99% bacteria => 99% bacteria, 0.5% hyena, 0.3% cow, 0.15% pig, 0.04% mammoth, 0.01% Human
# 0.5 -> 95912 0.3% -> 9913, 0.15% -> 9823, 0.04 -> 9785 0.01,  -> 9606

# tabl['A'].loc[30][tabl['A'].loc[30]>=0.0003].idxmin()

def get_sequence_with_substitution(sequence):
    """Generate a sequence with mutations relative to the given mutation probability matrix
    
    Parameters
    ----------
    sequence : Bio.Seq.MutableSeq
    
    Returns
    -------
    Bio.Seq.Seq
    """
    read_length = len(sequence)
    # create a set of draws
    choices = np.random.random(read_length)
    # create a profile sequence with as many profile 30 relative to read_length
    positions = np.concatenate([np.arange(15), np.arange(29, 14, -1)])
    positions = np.insert(positions, 15, np.full(read_length - 30, 30))
    newSeq = list(sequence)
    # iterate over the sequence, choices and positions
    for idx, (nc, (choice, pos)) in enumerate(zip(sequence, zip(choices, positions))):
        # we draw a probability which allows us to mutate
        if choice > same_nuc[pos][nc]:
            t = tabl[nc].loc[pos]
            try:
                newSeq[idx] = t[choice < t].idxmin()
            except:  # rounding error when over choice 0.999999
                newSeq[idx] = t.idxmax()
    return Seq(''.join(newSeq))

# specific code to work with mappability track...


def _get_map_pos(n_samples, map_file="/tmp/fred/map_track.bed.gz"):
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


def get_chrom(refseq):
    """Get the chromosome from a RefSeq entry.
    
    Parameters
    ----------
    refseq : Bio.Seq.Seq
        A RefSeq entry
    
    Returns
    -------
    str
        The chromosome extracted from the sequence id
    """
    chrom = refseq.id  # in case of scaffold genome
    if "|" in chrom:
        chrom = chrom.split("|")[0]
    return chrom

# sort the list of requences according to chromosome


def sort_recs(recs, chromosome_list=[], split_char="|"):
    """Sort a list of Bio.SeqRecord according to the chromosome, position, length.
    
    Parameters
    ----------
    recs : list of Bio.SeqRecord
    chromosome_list : list of str
        A list specifying the order of the chromosomes (as python doesn't really natural sort lists)
    split_char : str
        The delimiter used for splitting the record.id
    
    Returns
    -------
    list of Bio.SeqRecord
        Records properly sorted
    """
    def sort_func(item):
        """we sort according to chromosomes
         our header is >16|kraken:taxid|9906|69694935_57 ...
         use split to access the first element
         last element is pos_length
        """
        chrom = item.id.split(split_char)[0]
        pos, l = item.id.split(split_char)[-1].split("_")
        if chromosome_list:
            return (chromosome_list.index(chrom), int(pos), int(l))
        else:
            return (int(pos), int(l))
    return sorted(recs, key=sort_func)


def simulate_deamination(sequence, nbases=3):
    """Add C -> T mutations to 5' and 3' ends
    
    Parameters
    ----------
    sequence : Bio.Seq.MutableSeq
    nbases : int
        Number of bases to mutate for both ends
        
    Returns
    -------
    Bio.Seq.Seq
        The deaminated sequence
    """
    while "C" in sequence[:nbases]:  # 5' C to T
        sequence[sequence.index("C")] = "T"
    while "C" in sequence[-nbases:]:  # 3' C to T
        sequence[-nbases + sequence[-nbases:].index("C")] = "T"
    return sequence.toseq()


def mutate_unif(sequence, unif):
    """Add mutations uniformely throughout the sequence
    
    Parameters
    ----------
    sequence: Bio.Seq.Seq
    unif: float
        The mutaion rate

    Returns
    -------
    Bio.Seq.Seq
        A new mutated sequence
    """
    same_nuc = set("ACGT")
    read_length = len(sequence)
    # we will mutate a base only if our sample is < unif
    choices, = np.where(np.random.random(read_length) < unif)
    newSeq = list(sequence)
    for idx in choices:
        # we call tuple in order to do random.choice (and remove the current nc from the possibilities)
        newSeq[idx] = random.choice(tuple(same_nuc.difference(sequence[idx])))
    return Seq(''.join(newSeq))


# the idea is to generate a set of coordinates/size pairs and only take those substrings:
# sample 0:length_genome (N = #chunks desired)
# sample N sizes (N = chunks desired)
def chunk_fast(record, n_samples, vcf_in=None, chrom=None, specie="|", individual=0, unif=False, len_distrib=False, deaminate=0, minlength=35, maxlength=100, nthreads=4):
    """Split the sequence into chunks and add mutations (in parallel)
    
    The function will generate a list of positions to extract reads from as well as a list of read length (following the given distribution), it will then create a list of reads which will be mutated according to the mutation rule given (unif, VCF, or matrix)
    
    Parameters
    ----------
    record: Bio.SeqRecord
    n_samples: int
        How many reads should be sampled
    vcf_in: pysam.VariantFile, optional
        File providing mutations informations (VCF or Substitution matrix)
    chrom: int, optional
    specie: str, optional
    individual: int, optional
    unif: float, optional
        The mutaion rate, will be added on top of VCF/substitution Matrix 
    len_distrib: str, optional
    deaminate: int, optional
    minlength: int, optional
    maxlength: int, optional
    nthreads: int, optional
    
    Returns
    -------
    list of Bio.Seq.Seq
    """
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
        # length = random.choices(range(minlength, maxlength), k = n_samples) #only from python3.6
        # onwards
        length = np.random.choice(
            np.arange(minlength, maxlength+1), n_samples)
    samples = []
    for pos, l in zip(positions, length):
        rec = record[pos:pos + l]
        while not set(rec).issubset('ACGT'):  # control for unknown bases
            # print ("Unknown base (N) in read at {} {},
            # resampling...".format(pos,l), file=sys.stderr) # debug purpose
            pos = random.choice(range(0, len(record) - minlength))
            if len_distrib:
                l = np.random.choice(
                    np.arange(minlength, maxlength+1), 1, p=p)[0]
            else:
                l = random.choice(range(minlength, maxlength+1))
            rec = record[pos:pos + l]

        # create a new header which includes the specie read pos, read length
        samples += [SeqRecord.SeqRecord(rec.seq, id="{}{}{}_{}".format(rec.id, specie, pos, l),
                                        description=" ".join(rec.description.split(' ')[1:]))]
    if vcf_in:  # insert variation from VCF
        res = []
        if isinstance(vcf_in, pysam.VariantFile):  # we use a VCF
            for sample in samples:
                sequence = sample.seq.tomutable()
                l = len(sequence)
                for vcf_rec in vcf_in.fetch(chrom, start=pos, end=pos + l):
                    # Insert alt only if GT != 0/0, use the first individual by
                    # default
                    if any(vcf_rec.samples[0]['GT']):
                        if len(vcf_rec.alts[0]) == 1:  # only SNPs
                            sequence[vcf_rec.start - pos] = vcf_rec.alts[0]
                if deaminate or unif:  # save the seqs for later deam/unif
                    res.append(sequence)
                else:
                    sample.seq = sequence.toseq()
                # print(record[vcf_rec.start-3:vcf_rec.start+4].seq,vcf_rec.ref,vcf_rec.alts,vcf_rec.start,
                # pos, l) # debug purpose
        else:  # we use a variation probability matrix
            # run on multithreaded code as it is the bottleneck (each base is possibly mutated)
            with multiprocessing.Pool(processes=nthreads) as pool:
                res = pool.map(get_sequence_with_substitution, [
                               s.seq.tomutable() for s in samples])
        if len(res):  # run multi-threaded code
            with multiprocessing.Pool(processes=nthreads) as pool:
                if deaminate:
                    deam_func = partial(
                        simulate_deamination, deaminate=deaminate)
                    res = pool.map(deam_func, res)
                if unif:
                    unif_func = partial(mutate_unif, unif=unif)
                    res = pool.map(unif_func, res)
            for s, seq in zip(samples, res):
                s.seq = seq
    return samples


def read_len_distrib(filename, minlen=35, maxlen=100):
    """Reads a tsv file with counts of read per length
    
    Parameters
    ----------
    filename : str
        tsv file containing the read length count
    minlen : int
        minimum length for the distribution
    maxlen : int
        maximum length for the distribution
        
    Returns
    -------
    list of float
        The computed proportion of each read length
    """
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
    """Reads a fasta file and computes the number of sequences to extract per chromosome (relative to chromosome length)
    
    Parameters
    ----------
    file_in : str
        FASTA file you want to estimate the distribution from.
    num_seq : int
        The total number of sequences you want to sample.
    n_chromosomes : int
        If you want to take only the first *n* chromosomes of your dataset.
    
    Returns
    -------
    list of int
        The number of sequences to extract per chromosome.
    """
    try:
        with pysam.FastaFile(file_in) as fa:
            print("The file contains", fa.nreferences,
                  "sequences", file=sys.stderr)
            if num_seq:
                if n_chromosomes:
                    print("Working only with the first",
                          n_chromosomes, file=sys.stderr)
                    full_size = sum(fa.lengths[:n_chromosomes])
                    reads_per_chrom = [
                        int(s / full_size * num_seq) + 1 for s in fa.lengths[:n_chromosomes]]
                else:
                    n_chromosomes = len(fa.lengths)
                    full_size = sum(fa.lengths)
                    reads_per_chrom = [
                        int(s / full_size * num_seq) + 1 for s in fa.lengths]
    except OSError:  # fasta is not indexed do a naive fallback
        print("Error, fasta file not indexed, doing naive sampling...",
              file=sys.stderr)
        if n_chromosomes:
            reads_per_chrom = [
                int(num_seq / n_chromosomes) + 1] * n_chromosomes
        else:
            print(
                "Naive sampling needs the option --chromosomes...Aborting...", file=sys.stderr)
            sys.exit(1)
    # our samplig tends to slightly over sample, readjust by removing reads on random chromosomes
    extra_samples = sum(reads_per_chrom) - num_seq
    for idx in [random.randint(0, n_chromosomes-1) for _ in range(extra_samples)]:
        # some scaffolds are so small that the don't even get a read to sample.
        while reads_per_chrom[idx] == 0:
            idx = random.randint(0, n_chromosomes-1)
        reads_per_chrom[idx] -= 1
    return reads_per_chrom


# dummy function used for opening a file non-gzipped
@contextlib.contextmanager
def _ret_file(f):
    yield f

def read_mutation_matrix(file_in):
    """Create a table from a mutation matrix provided by SNPad
    
    Returns
    -------
    tbl : Pandas.table
        A table containing the mutation probability for base at each position of a read.
    sn : list of dict
        The list contains the highest probability for which we can have a mutation on a given position.
        Used to optimize the computation when we draw for a mutation.
    """
    t = pd.read_table(file_in, sep="\t",
                      header=None, names=["profile", "from", "to", "prob"])
    tabl = pd.pivot_table(t, values='prob', index=[
                          'profile', 'to'], columns=['from'], aggfunc=np.sum)
    # contains, for each nc, the highest probability for which we can have
    # a mutation
    # e.g.
    # 'A': 0.00797800000000004, -> ~0.79% chances to have a mutation, any sampled number above this would result in keeping the nc
    # this allows us to optimize the algorithm by not using the lookup
    # table and directly keep the nc.
    sn = [{nc: prob for nc, prob in zip(tabl.columns,
                                                  [tabl[x][pos][x] for x in tabl])} for pos in range(31)]
    return tbl, sn

def _main():
    random.seed()
    # Process command line
    parser = argparse.ArgumentParser(description='Split a genome into chuncks')
    parser.add_argument('--num_seq', default=0, type=int,
                        help='Number of sequences to output')
    parser.add_argument(
        '--outfile',  help='Fasta file where to extract sequence')
    parser.add_argument('--specie',  help='Specie Taxa')
    parser.add_argument(
        'file_in', metavar="genome.fa[.gz], [bgzip] and indexed")
    parser.add_argument('--chromosomes', type=int,
                        help='How many chromosomes are expected', default=0)
    parser.add_argument('--minlen', type=int,
                        help='Minimum read length (lengths are sampled uniformly)', default=35)
    parser.add_argument('--maxlen', type=int,
                        help='Maximum read length (lengths are sampled uniformly)', default=100)
    parser.add_argument('--threads', type=int, default=multiprocessing.cpu_count(),
                        help="Specify number of threads to use")
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
    sorting = parser.add_mutually_exclusive_group()
    sorting.add_argument('--sorted', action='store_true',
                         help='Natural sort the reads by chromosome name, position and length')
    sorting.add_argument('--shuffled', action='store_true',
                         help='Output shuffled reads')
    args = parser.parse_args()

    print("Loading dataset...", end="", file=sys.stderr)

    if args.substitution_file:
        global tabl
        global same_nuc
        tabl, same_nuc = read_mutation_matrix(args.substitution_file)
        for nc in 'ACGT':
            for pos in range(31):
                # do cumsum on sorted values, assuming that the highest proba is the nc itself: C->C, A->A, etc
                tabl[nc][pos] = [tabl[nc][pos].sort_values(ascending=False).cumsum()[
                    to] for to in tabl.columns]
    num_reads_to_sample = estimate_read_distribution(
        args.file_in, args.num_seq, args.chromosomes)
    specie = "|"
    if args.specie:
        specie += args.specie.replace(' ', "_")+"|"
    with gzip.open(args.file_in, "rt") if args.file_in.endswith("gz") else _ret_file(args.file_in) as file_in, open(args.outfile, 'w') if args.outfile else sys.stdout as file_out:
        all_chunks = []
        vcf_in = args.substitution_file  # will be none if we use a VCF file
        chromosomes = []
        if args.vcf:
            vcf_in = pysam.VariantFile(args.vcf)
            # subset the VCF to 1 individual only
            vcf_in.subset_samples([vcf_in.header.samples[args.individual]])
            # our chromosome can contain up to 3 characters (most likely 2 max:
            # e.g. chromosome 21,)
        p = re.compile("chromosome \w{1,3},", re.IGNORECASE)
        for num_record, record in enumerate(SeqIO.parse(file_in, "fasta")):
            # we want reads only to assigned chromosomes
            if args.chromosomes:
                if num_record > args.chromosomes - 1:  # num_record is 0-based
                    # assume the genome is sorted, with chromosomes first...
                    break
                print("Parsing chromosome", num_record +
                      1, end="\r", file=sys.stderr)
            if num_reads_to_sample[num_record] == 0:  # skip this sequence
                continue
            chrom = p.search(record.description)
            if not chrom:
                chrom = get_chrom(record)  # in case of scaffold genome
            else:
                chrom = chrom.group()[len('chromosome '):-1]
            if args.vcf or args.substitution_file:
                # we write chromosomes in the same order as the given input.
                if not args.shuffled and not args.sorted:
                    all_chunks = sort_recs(chunk_fast(record, num_reads_to_sample[
                        num_record], vcf_in, chrom, specie, unif=args.unif, deaminate=args.deaminate, len_distrib=args.length,
                        minlength=args.minlen, maxlength=args.maxlen, nthreads=args.threads))
                    SeqIO.write(all_chunks, file_out, "fasta")
                else:
                    chromosomes.append(chrom)
                    all_chunks += chunk_fast(record, num_reads_to_sample[
                        num_record], vcf_in, chrom, specie, unif=args.unif, deaminate=args.deaminate, len_distrib=args.length, minlength=args.minlen, maxlength=args.maxlen, nthreads=args.threads)

            if num_record % 100 == 99:  # show progress
                print(num_record + 1, "sequences parsed...",
                      end="\r", file=sys.stderr)
        if args.vcf:
            vcf_in.close()
        print("Done\nWritting down records...", end="", file=sys.stderr)
        if args.shuffled:
            random.shuffle(all_chunks)
        if args.sorted:  # we will do a natural sort on the chromosomes
            chromosomes.sort(key=lambda key: [int(text) if text.isdigit(
            ) else text for text in re.split('([0-9]+)', key)])
            all_chunks = sort_recs(all_chunks, chromosome_list=chromosomes)
        if args.shuffled or args.sorted:
            SeqIO.write(all_chunks, file_out, "fasta")
    print("Done", file=sys.stderr)


if __name__ == "__main__":
    _main()
