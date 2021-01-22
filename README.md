# metagen utils
a set of tools to work with ancient metagenomics data
## installation
pip install --user metagen_utils-VERSION.tar.gz
- or -
unpack
cd /WHERE_UNPACKED

python -m metagen_utils-VERSION.chunk_genome -h
## generate the doc
unpack
cd metagen_utils-VERSION
python setup.py build_sphinx

## Data
files present in data/ are are empirical values which enable `chunk_genome` to create simulated reads close to what one would expect in reality regarding degradation of DNA and read length distribution.

kraken_report needs trimmed data from the NCBI taxonomy ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz put in you kraken taxonomy folder
On UNIX:
grep 'scientific name' names.dmp |cut -d'|' -f 1,2 |gzip -c > /KRAKEN_DB/taxonomy/names_trimmed.dmp.gz
cut -d '|' -f 1,2,3 nodes.dmp|gzip -c > /KRAKEN_DB/taxonomy/nodes_trimmed.dmp.gz
