import os
#snakemake -p -j 4 -s AMPS.Snakefile --config threads=4 infile="~/Projects//Mash/simulationA2_99bacteria.fasta"

filtr = config["filter"] if "filter" in config else 'def_anc'
infile = config["infile"]
database = config["db"] if "db" in config else "/home/frederic_romagne/Projects/malt-git/refseq_kraken_malt"
specie_list = config["specie"] if "specie" in config else 'Homo,Homo_sapiens_ssp._Denisova'
mpi = config["mpi"] if "mpi" in config else '85.0'
sup = config["sup"] if "sup" in config else '1'
mq = config["mq"] if "mq" in config else '100'
top = config["top"] if "top" in config else '1'
nthread = config["threads"] if "threads" in config else '1'
malt_run = 'malt-git/malt-run'
java='/home/frederic_romagne/jdk-9.0.4/bin/java'

snakepath = os.path.dirname(workflow.snakefile)

#resources is a wildcard
reso = snakepath+'/resources/'

def filter_output():
    if filtr == "def_anc":
        return ["maltExtract/default/RunSummary.txt",
        "maltExtract/default/TotalCount.txt",
        "maltExtract/ancient/RunSummary.txt",
        "maltExtract/ancient/TotalCount.txt"]
    else:
        return ["maltExtract/"+filtr+"/RunSummary.txt",
        "maltExtract/"+filtr+"/TotalCount.txt"]

rule all:
   #java -cp /home/frederic_romagne/NetBeansProjects/AMPS/dist/AMPS.jar:\
   #/home/frederic_romagne/NetBeansProjects/RMASifterResurgencePrime/dist/RMASifterResurgencePrime.jar:\
   #/home/frederic_romagne/Projects/jloda-git/antbuild/class/:\
   #'/home/frederic_romagne/Projects/malt-git/antbuild/MALT.jar':\
   #'/home/frederic_romagne/Projects/rmasifter-the-rmaing/RMASifterResurgencePrime/dependencies/commons-cli-1.3.1/commons-cli-1.3.1.jar' Main --configFile amps.conf \
   #-i ../Mash/simulationA2_99bacteria.fasta -o toto
    input:
        "malt/",
        filter_output(),
        'maltExtract/pdf_candidate_profiles/'

rule Malt:
    ### cp = /home/frederic_romagne/Projects/malt-git/jars/.install4j/i4jruntime.jar:\
    #   /home/frederic_romagne/Projects/malt-git/class:\
    #   /home/frederic_romagne/Projects/malt-git/jars/data.jar
    #
    #CMD
    #/home/frederic_romagne/Projects/malt-git/malt-run \
    #-d /home/frederic_romagne/Projects/malt-git/refseq_kraken_malt \
    #-i /home/frederic_romagne/Projects/Mash/simulationA2_99bacteria.fasta \
    #-o /home/frederic_romagne/Projects/amps/toto/malt \
    #-m BlastN -at SemiGlobal --memoryMode load -t 1 -sup 1 -mq 100 -top 1 -mpi 85.0 -v ""
    input:
        db = database,
        seq = infile
    output:
        "malt/"
    shell:
        "{malt_run} \
        -d {input.db} \
        -i {input.seq} \
        -o malt \
        -m BlastN -at SemiGlobal --memoryMode load -t {nthread} -sup {sup} -mq {mq} -top {top} -mpi {mpi}"


rule MaltExtract:
    #java -cp /home/frederic_romagne/NetBeansProjects/RMASifterResurgencePrime/dist/RMASifterResurgencePrime.jar:\
    #/home/frederic_romagne/Projects/jloda-git/antbuild/class/:\
    #'/home/frederic_romagne/Projects/rmasifter-the-rmaing/RMASifterResurgencePrime/dependencies/commons-math3-3.6.1.jar':\
    #'/home/frederic_romagne/Projects/rmasifter-the-rmaing/RMASifterResurgencePrime/dependencies/forester_1041.jar':\
    #'/home/frederic_romagne/Projects/malt-git/antbuild/MALT.jar':\
    #'/home/frederic_romagne/Projects/rmasifter-the-rmaing/RMASifterResurgencePrime/dependencies/commons-cli-1.3.1/commons-cli-1.3.1.jar'\
    # RMAExtractor.RMAExtractor -i /home/frederic_romagne/Projects/amps/toto/malt/ -o /home/frederic_romagne/Projects/amps/toto/maltExtract/ -p 1 -t Homo_sapiens -f def_anc --resources /home/frederic_romagne/Projects/rmasifter-the-rmaing/RMASifterResurgencePrime/resources/ --top 0.01
    input:
        "malt/"
    output:
        filter_output()
    shell:
        "{java} -cp {snakepath}/RMASifterResurgencePrime.jar:{snakepath}/jloda-git/class/:{snakepath}/'dependencies/commons-math3-3.6.1.jar':{snakepath}/'dependencies/forester_1041.jar':{snakepath}/'malt-git/jars/MALT.jar':{snakepath}/'dependencies/commons-cli-1.3.1/commons-cli-1.3.1.jar' \
     RMAExtractor.RMAExtractor -i {input} -o maltExtract \
     -p {nthread} -t {specie_list} -f {filtr} --top 0.01 \
     --resources {reso} -v"
     #/home/frederic_romagne/Projects/rmasifter-the-rmaing/RMASifterResurgencePrime/resources/

rule PostProc:
    input:
        'maltExtract/'
    output:
        'maltExtract/pdf_candidate_profiles/'
    shell:
        "postproc/postprocessing.AMPS.r -m {filtr} -r {input} -t {threads} -n {specie_list}"
