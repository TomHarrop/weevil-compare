#!/usr/bin/env python3

import csv
import pathlib as pathlib2
import tempfile
import multiprocessing

###########
# GLOBALS #
###########

fn_file = 'data/assembly_filenames.txt'
assembly_path = 'data/assemblies'
l1r1 = 'data/reads/CCU1EANXX-3897-01-21-1_S123_L005_R1_001.fastq.gz'
l1r2 = 'data/reads/CCU1EANXX-3897-01-21-1_S123_L005_R2_001.fastq.gz'
l2r1 = 'data/reads/CCU1EANXX-3897-02-21-1_S124_L005_R1_001.fastq.gz'
l2r2 = 'data/reads/CCU1EANXX-3897-02-21-1_S124_L005_R2_001.fastq.gz'

busco_container = 'shub://TomHarrop/singularity-containers:busco_3.0.2'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
star_container = 'shub://TomHarrop/singularity-containers:star_2.7.0c'
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
samtools_container = 'shub://TomHarrop/align-utils:samtools_1.10'
r_container = 'shub://TomHarrop/singularity-containers:r_3.5.2'
minimap2 = 'shub://TomHarrop/singularity-containers:minimap2_2.11r797'


#############
# FUNCTIONS #
#############


def assembly_catalog_resolver(wildcards):
    if wildcards.name in spec_to_file:
        return({'fasta': spec_to_file[wildcards.name]})
    else:
        raise ValueError('missing {} in catalog'.format(wildcards.name))


def genomic_read_resolver(wildcards):
    if wildcards.mode == 'sr':
        return('data/illumina_reads.fq.gz')
    elif wildcards.mode == 'map-ont':
        return('data/nanopore_reads.fq.gz')
    else:
        raise ValueError(f'genomic_read_resolver can\'t resolve {wildcards.mode}')


def resolve_path(x):
    return str(pathlib2.Path(x).resolve())


########
# MAIN #
########


spec_to_file = {}
spec_to_commonname = {}
with open(fn_file, 'rt') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        spec_to_file[row[1]] = str(
            pathlib2.Path(assembly_path, row[0]))
        spec_to_commonname[row[1]] = row[3]


rnaseq_samples = glob_wildcards(
    'data/reads/CCU1EANXX-3897-{s1}-21-1_S{s2}_L005_R2_001.fastq.gz')
sample_names = expand('CCU1EANXX-3897-{s1}-21-1_S{s2}',
                      zip,
                      s1=rnaseq_samples.s1,
                      s2=rnaseq_samples.s2)

#########
# RULES #
#########

rule target:
    input:
        expand('output/010_busco/run_{name}/full_table_{name}.tsv',
               name=list(spec_to_file.keys())),
        'output/020_stats/stats.txt',
        'output/030_map/star_results.csv',
        expand('output/040_gbs_map/{name}.flagstat',
               name=list(spec_to_file.keys())),
        expand('output/050_map-genomic-reads/{mode}/{name}.flagstat',
               name=list(spec_to_file.keys()),
               mode=['map-ont', 'sr'])

# genomic read mapping
rule map_genomic_reads_stats:
    input:
        'output/050_map-genomic-reads/{mode}/{name}_sorted.bam'
    output:
        'output/050_map-genomic-reads/{mode}/{name}.flagstat'
    singularity:
        samtools_container
    priority:
        10
    shell:
        'samtools flagstat -O tsv {input} > {output}'


rule minimap:
    input:
        reads = genomic_read_resolver,
        ref = 'output/050_map-genomic-reads/minimap2-index_{mode}/{name}/ref.mmi'
    output:
        'output/050_map-genomic-reads/{mode}/{name}.sam'
    log:
        'output/logs/050_map-genomic-reads/{mode}_{name}.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        minimap2
    shell:
        'minimap2 '
        '-ax {wildcards.mode} '
        '-t {threads} '
        '{input.ref} '
        '{input.reads} '
        '> {output} '
        '2> {log}'


rule minimap2_index:
    input:
        unpack(assembly_catalog_resolver)
    output:
        'output/050_map-genomic-reads/minimap2-index_{mode}/{name}/ref.mmi'
    log:
        'output/logs/050_map-genomic-reads/{mode}_{name}_index.log'
    wildcard_constraints:
        name = '\w+',
        mode = '(\w|-)+'
    singularity:
        minimap2
    shell:
        'minimap2 '
        '-x {wildcards.mode} '
        '-d {output} '
        '{input.fasta} '
        '2> {log}'

# catalog mapping
rule map_stacks_catalog_stats:
    input:
        'output/040_gbs_map/{name}.sam'
    output:
        'output/040_gbs_map/{name}.flagstat'
    singularity:
        samtools_container
    priority:
        10
    shell:
        'samtools flagstat -O tsv {input} > {output}'

rule map_stacks_catalog:
    input:
        fa = 'data/catalog.fa.gz',
        index = expand('output/040_gbs_map/bwa_index_{{name}}/ref.fasta.{suffix}',
                       suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        'output/040_gbs_map/{name}.sam'
    params:
        prefix = 'output/040_gbs_map/bwa_index_{name}/ref.fasta',
    threads:
        min(multiprocessing.cpu_count(), 48)
    log:
        'output/logs/map_stacks_catalog_{name}.log'
    singularity:
        bwa_container
    shell:
        'bwa mem '
        '-t {threads} '
        '{params.prefix} '
        '{input.fa} '
        '> {output} '
        '2> {log}'

rule bwa_index:
    input:
        unpack(assembly_catalog_resolver)
    output:
        expand('output/040_gbs_map/bwa_index_{{name}}/ref.fasta.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    params:
        prefix = 'output/040_gbs_map/bwa_index_{name}/ref.fasta'
    threads:
        1
    log:
        'output/logs/bwa_index_{name}.log'
    singularity:
        bwa_container
    shell:
        'bwa index '
        '-p {params.prefix} '
        '{input.fasta} '
        '2> {log}'


# RNAseq read mapping
rule parse_star_results:
    input:
        star_files = expand('output/030_map/{s}_{name}/{s}.Log.final.out',
                            s=sample_names,
                            name=list(spec_to_file.keys())),
        sample_key = 'data/full_sample_key.csv',
        assembly_filenames = 'data/assembly_filenames.txt'
    output:
        'output/030_map/star_results.csv'
    params:
        star_dir = 'output/030_map'
    log:
        'output/logs/030_map/parse_star_results.log'
    singularity:
        r_container
    script:
        'src/parse_star_results.R'


rule map:
    input:
        r1 = 'data/reads/{s}_L005_R1_001.fastq.gz',
        r2 = 'data/reads/{s}_L005_R2_001.fastq.gz',
        genome = 'output/030_map/star-index_{name}/SA'
    output:
        ('output/030_map/{s}_{name}/'
         '{s}.Log.final.out')
    params:
        genome_dir = 'output/030_map/star-index_{name}',
        prefix = 'output/030_map/{s}_{name}/{s}.'
    threads:
        min(multiprocessing.cpu_count(), 48)
    log:
        'output/logs/030_map/{s}_{name}.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--readFilesCommand zcat '
        '--outSAMtype BAM Unsorted '
        '--outBAMcompression 10 '
        '--readFilesIn {input.r1} {input.r2} '
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'

rule generate_genome:
    input:
        unpack(assembly_catalog_resolver)
    output:
        'output/030_map/star-index_{name}/SA',
        'output/030_map/star-index_{name}/chrNameLength.txt'
    params:
        outdir = 'output/030_map/star-index_{name}'
    threads:
        multiprocessing.cpu_count()
    log:
        'output/logs/030_map/generate_genome_{name}.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--limitGenomeGenerateRAM 500000000000 ' # 0.5 TB
        '--genomeDir {params.outdir} '
        '--genomeFastaFiles {input.fasta} '
        '--outFileNamePrefix {params.outdir}/ '
        '&> {log}'


# stats
rule stats:
    input:
        fasta = list(spec_to_file.values())
    output:
        stats = 'output/020_stats/stats.txt'
    params:
        fasta = lambda wildcards, input: ','.join(input.fasta)
    log:
        'output/logs/020_stats/statswrapper.log'
    threads:
        min(len(list(spec_to_file.values())),
            multiprocessing.cpu_count())
    singularity:
        bbduk_container
    shell:
        'statswrapper.sh '
        'in={params.fasta} '
        'minscaf=1000 '
        'format=3 '
        'threads={threads} '
        '> {output.stats} '
        '2> {log}'


rule busco:
    input:
        unpack(assembly_catalog_resolver),
        lineage = 'data/endopterygota_odb9'
    output:
        'output/010_busco/run_{name}/full_table_{name}.tsv'
    log:
        resolve_path('output/logs/010_busco/busco_{name}.log')
    params:
        wd = 'output/010_busco',
        fasta = lambda wildcards, input: resolve_path(input.fasta),
        lineage = lambda wildcards, input: resolve_path(input.lineage),
        tmpdir = tempfile.mkdtemp()
    threads:
        multiprocessing.cpu_count()
    singularity:
        busco_container
    shell:
        'cd {params.wd} || exit 1 ; '
        'run_BUSCO.py '
        '--force '
        '--tmp_path {params.tmpdir} '
        '--in {params.fasta} '
        '--out {wildcards.name} '
        '--lineage {params.lineage} '
        '--cpu {threads} '
        '--species tribolium2012 '
        '--mode genome '
        '&> {log}'


rule samtools_sort:
    input:
        '{path}/{file}.sam'
    output:
        bam = '{path}/{file}_sorted.bam',
        bai = '{path}/{file}_sorted.bam.bai'
    log:
        'output/logs/sort/{path}_{file}.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        samtools_container
    shell:
        'samtools sort '
        '-o {output.bam} '
        '-O BAM '
        '-@ {threads} '
        '{input} '
        '2> {log} '
        '; '
        'samtools index {output.bam} 2>> {log}'

