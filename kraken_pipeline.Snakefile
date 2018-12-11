#~/.local/bin/snakemake --jobs 200 -s /mnt/sediments/fred/metagen.p3.kraken.SnakeFile --cluster-config /mnt/sediments/fred/sediment_config.json --cluster "qsub -pe smp {threads} -l class=cow" --local-cores 48 --latency-wait 30 --restart-times 3

#grep '\sancient' out/kraken/*/substitution/BinomTest_CtoT.txt|cut -d. -f 2|cut -d/ -f5|sort|uniq

#specie_id = {'Primates': [('9606', 'Homo_sapiens')], 'Bovidae': [('9606', 'Homo_sapiens')]}

#import random
#def get_specie_id_name (family_name):
#    return random.choice(specie_id[family_name])
# remove_duplicate query_length = 35

rm_dup_query_len = config["rm_dup_query_len"] \
if "rm_dup_query_len" in config else 35
num_min_dup = int(config["num_min_dup"]) if "num_min_dup" in config else 2

    
localrules: all, demultiplex, remove_duplicate, kraken 
#ruleorder:  kraken_report > bwa_map



rule all:
    input:
        dynamic("{libID}.kraken.done"),#,zip,libID=['SP3593','SP3595'],'pseudouniq/pseudouniq_stats.txt' family_name=['Primates','Bovidae']))
        #'pseudouniq/pseudouniq_stats.txt'
#        dynamic("out/kraken/{family_name}/aligned/{libID}.bam"),
#        dynamic('split/{libID}.bam'),
#        dynamic('{libID}.kraken.report'),
#        "most_clades.{libID}.txt",
#        'pseudouniq/{libID}.noPCRdups.fq'
        #"out/kraken/{family_name}/aligned/{libID}.{specie_name}.bam"


rule demultiplex:
    input:
        bam = config["bamfile"],
        byfile = config["byfile"]
    output:
        dynamic('split/{libID}.bam')
    threads: 1
    shell:
        'python3 ~frederic_romagne/matthias_python_reimplementations/split_bam.py --byfile {input.byfile} --outdir split --minscore 10 --max 0 {input.bam} > split/splittingstats.txt'

rule remove_duplicate:
    input:        
        bam="split/{libID}.bam"
    output:
        bam="pseudouniq/{libID}.noPCRdups.bam",
        stats="pseudouniq/{libID}.pseudouniq_stats.txt"
    threads: 1
    script:
        "bam-rmdup.py"

rule merge_pseudouniq_stats:
    input:
        rules.remove_duplicate.output.stats
    output:
        'pseudouniq/pseudouniq_stats.txt'
    shell:
        'echo -e "#file\tall_seqs\tseqs>={rm_dup_query_len}\t\
            avrg_times_seen_L>={rm_dup_query_len}\tfinal_noPCRdups_seqs\t\
            avrg_times_seen_final" > {output}\
           cat {input} | sort | uniq >> {output}'

rule bam2fastq:
    input:
        #dynamic('pseudouniq/{libID}.noPCRdups.bam')
        'pseudouniq/{libID}.noPCRdups.bam'
    output:
        'pseudouniq/{libID}.noPCRdups.fq'
    threads: 1
    shell:
        'bam2fastx -q -Q -A -o {output} {input}'
        
rule kraken:
    input:
        'pseudouniq/{libID}.noPCRdups.fq'
    output:
        '{libID}.kraken'
    threads: 15
    shell: # run kraken on 15 threads in order to parallelize it.
        '/home/frederic_romagne/kraken/install/kraken --threads {threads} --db /mnt/ramdisk/refseqReleaseKraken --only-classified-output --output {output} {input} 2>{output}.stats'

rule most_clade:
    input:
        '{libID}.kraken'      
    output:
        #most_clades=temp("most_clades.{libID}.txt"),
        report='{libID}.kraken.report'
    group: 'kraken'
    threads: 1
    shell:
        'python3 ~frederic_romagne/sediment_shotgun/kraken_report.py \
        --db /mnt/sequencedb/Refseq/refseqReleaseKraken {input} > {output.report}'
        
rule kraken_report:
    input:
        bamfile='pseudouniq/{libID}.noPCRdups.bam',
        kraken='{libID}.kraken',
        most_clades='{libID}.kraken.report'
    threads: 1
    group: 'kraken'
    params:
        outdir='out/kraken'      
    output:
        report=temp('{libID}.kraken.done')
    shell:        
        # use all the families where species are available from RefSeq except Hominidae (9604) 
        # as we extract Primates automatically and and Orycteropodidae (9816) as we cannot map to genome (too big!)
        'python3 ~frederic_romagne/sediment_shotgun/kraken_report.py --db /mnt/sequencedb/Refseq/refseqReleaseKraken --clade $(grep "\sF\s" SP6024.kraken.report| sort -rg -k2,2 | grep -f /mnt/sediments/fred/refseq_family_list.txt | head -n5 | cut -f5 | tr "\n" ,)9443 --extractFile {input.bamfile} --outdir {params.outdir} {input.kraken} > {output.report}'


onsuccess:
    shell('echo -e "#file\tall_seqs\tseqs>={rm_dup_query_len}\t\
            avrg_times_seen_L>={rm_dup_query_len}\tfinal_noPCRdups_seqs\t\
            avrg_times_seen_final" > pseudouniq/pseudouniq_stats.txt;\
           cat pseudouniq/*.pseudouniq_stats.txt | sort | uniq >> pseudouniq/pseudouniq_stats.txt')
        
#def get_families_from_sample(wildcards):
#    with open('most_clades.'+wildcards.libID+'.txt', 'r') as most_clades:    
#        return ['out/kraken/{}.{}.bam'.format(wildcards.libID, family_name[line.strip()]) for line in most_clades.readlines()]
        
#rule bwa_map:
#    input:
#        #report='{libID}.kraken.report',
#        #assigned=get_families_from_sample
#        assigned=dynamic('out/kraken/{libID}.{family_name}.bam')
#        #assigned=dynamic('out/kraken/{libID}.assigned.{family_name}.bam')
#        #aligned=rules.kraken_report.output
#    threads: 8
#    params:
#       specie_id=lambda wildcards: random.choice(specie_id[wildcards.family_name])
#    output:
#        bamfile='out/kraken/{family_name}/aligned/{libID}.bam'
#        #bamfile='out/kraken/{family_name}/aligned/{libID}.{specie_name}.bam'
#    shell:# generate specie_id and specie_name from rule kraken_report
#        '/home/public/usr/bin/bwa bam2bam -t {threads} -g /mnt/sequencedb/Refseq/tmp_ncbi/vertebrate_mammalian/{params.specie_id}.fa.gz -n 0.01 -o 2 -l 16500 --temp-dir /dev/shm --only-aligned {input.assigned} | samtools sort -@30 -o {output.bamfile}'
#        


