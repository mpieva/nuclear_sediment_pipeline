# metagen utils
This consists of 2 tools for analysis of metagenomic DNA from sediments.
- **kraken_report** is a rewrite in python of the original **kraken-report.pl** perl script.
  - You can now fetch back reads assigned by **kraken**, for user-defined clades
  - Improved speed

- **chunk_genome** provides functions used to generate reads derived from a set of sequences (i.e. full genome FASTA file).
  - You can specify a mutation rate, a VCF file in order to add known SNPs, or a mutation matrix obtained from snpAD
  > Kay Prüfer, snpAD: an ancient DNA genotype caller, Bioinformatics, Volume 34, Issue 24, 15 December 2018, Pages 4165–4171, https://doi.org/10.1093/bioinformatics/bty507.
  - Default read length distribution is uniform, but you can provide your own (see [Data](#data)).
  - Work is processed in parallel in order to decrease computation time.

## Usage
### Examples
- Report the clade assignment from the output of `kraken` and extract the sequences for clade 9606

  `python3 -m metagen_utils.kraken_report --db ${KRAKEN_DB_PATH} --clades 90606 --outdir results \`
  
  `--extractFile sourcefile.fa sourcefile.fa.kraken`
- Generate 30 reads with a specific deamination pattern only using the 1st chromosome of the source_genome

  `python3 metagen_utils.chunk_genome --subst ${DATA}/substitution_matrix.tsv  --num_seq 30 \`
  
  `--chrom 1 --outfile simulated_reads.fa source_genome.fa`

## Installation
`pip install --user metagen_utils-VERSION.tar.gz`

`python -m metagen_utils.chunk_genome -h`

\- or -

unpack

`cd /WHERE_UNPACKED`

`python -m metagen_utils-VERSION.chunk_genome -h`

### <a id="data"></a>Data
Files present in `data/` are are empirical values which enable `chunk_genome` to create simulated reads close to what one would expect in reality regarding degradation of DNA and read length distribution. Copy them for further use as they wont be installed automatically.

**kraken_report** needs trimmed data from the NCBI taxonomy ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz put in you kraken taxonomy folder

On UNIX:

`grep 'scientific name' names.dmp |cut -d'|' -f 1,2 |gzip -c > ${KRAKEN_DB_PATH}/taxonomy/names_trimmed.dmp.gz`

`cut -d '|' -f 1,2,3 nodes.dmp|gzip -c > ${KRAKEN_DB_PATH}/taxonomy/nodes_trimmed.dmp.gz`
## generate the doc using SPHINX
unpack

`cd metagen_utils-VERSION`

`python setup.py build_sphinx`
