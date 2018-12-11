from pathlib import Path
from argparse import ArgumentParser
from pysam import AlignmentFile
from Bio.Seq import Seq
from collections import defaultdict


def main(bamin, bamout, libID, statsout, rm_dup_query_len, num_min_dup):
    with open(statsout,'w') as stats, AlignmentFile(bamin, "rb", check_sq=False) as samfile:
            #using template=samfile to copy the header
        with AlignmentFile(bamout, "wb", template=samfile) as bamout:
            reads_to_write = []
            duplicate = defaultdict(int)
            length_passed = 0
            passed = 0
            for allreads, aln in enumerate(samfile.fetch(until_eof=True), 1):
                if aln.query_length >= rm_dup_query_len: #35
                    length_passed += 1
                    if aln.is_reverse: # if the read is reverse
                        # the sequence has to be revcomp'ed'
                        seq = Seq(aln.query_sequence).reverse_complement()#create a sequence object
                    else:
                        seq = aln.query_sequence
                    duplicate[seq] += 1
                    #keep only reads apearing at least twice
			        #twice by defaul, now can be assigned through command line
			        # we want the sequence to be printed only once it reaches the min
			        #number of duplicates
                    if duplicate[seq] == num_min_dup:
                        passed += 1
                        # TODO, take the best score read
                        reads_to_write.append(aln)
            for read in reads_to_write: # add the XP tag before writing
                if read.is_reverse:
                    seq = Seq(read.query_sequence).reverse_complement()
                else:
                    seq = read.query_sequence
                read.set_tag('XP', duplicate[seq], value_type='i')
                bamout.write(read)
        #write stats
        dups_average = [v for v in duplicate.values() if (v>=num_min_dup)]
        all_average = [v for v in duplicate.values()]
        #avoid division by 0
        if len(all_average):
            rounded_avrg_all = "{:.1f}".format(sum(all_average)/len(all_average))
        else:
            rounded_avrg_all = "0"
        if len(dups_average):
            rounded_avrg_final =" {:.1f}".format(sum(dups_average)/len(dups_average))
        else:
            rounded_avrg_final = "0"
        print (libID, allreads, length_passed, rounded_avrg_all,\
         passed, rounded_avrg_final, sep="\t", file=stats) #print out the stats
         

if __name__ == "__main__":
    try:# called from Snakemake
        rm_dup_query_len = snakemake.config["rm_dup_query_len"] \
        if "rm_dup_query_len" in snakemake.config else 35
        num_min_dup = int(snakemake.config["num_min_dup"]) if "num_min_dup" in snakemake.config else 2
        main(snakemake.input.bam, snakemake.output.bam, snakemake.wildcards.libID, snakemake.output.stats, rm_dup_query_len, num_min_dup)
    except NameError:
        # Process command line
        parser = ArgumentParser(description='Remove (perfect) duplicates, and XP tag')
        parser.add_argument(
            '--outdir',  help='Fasta file where to extract sequence', default="pseudouniq")
        parser.add_argument(
            'bamin', metavar="in.bam", help="BAM file to dedup")
        parser.add_argument('--minlen', type=int,
                            help='Minimum read length', default=35)
        parser.add_argument('--num_min_dup', type=int,
                            help='Minimum time we see this sequence before marking as duplicate', default=2)           
        args = parser.parse_args()
        bamin = Path(args.bamin)
        bamout = Path(args.outdir, bamin.stem+'.noPCRdups'+bamin.suffix)
        outstats = Path(args.outdir, bamin.stem+'.pseudouniq_stats.txt')
        main(args.bamin, bamout, bamin.stem, outstats, args.minlen, args.num_min_dup)
