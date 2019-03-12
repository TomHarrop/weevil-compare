#!/usr/bin/env python3

import csv
import pathlib2
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

#############
# FUNCTIONS #
#############


def assembly_catalog_resolver(wildcards):
    if wildcards.name in spec_to_file:
        return({'fasta': spec_to_file[wildcards.name]})
    else:
        raise ValueError('missing {} in catalog'.format(wildcards.name))


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


#########
# RULES #
#########

rule target:
    input:
        expand('output/010_busco/run_{name}/full_table_{name}.tsv',
               name=list(spec_to_file.keys())),
        'output/020_stats/stats.txt',
        expand('output/030_map/{s}/{s}.Log.final.out',
               s=['CCU1EANXX-3897-01-21-1_S123',
                  'CCU1EANXX-3897-02-21-1_S124'])

# RNAseq read mapping
rule map:
    input:
        r1 = 'data/reads/{s}_L005_R1_001.fastq.gz',
        r2 = 'data/reads/{s}_L005_R2_001.fastq.gz',
        genome = 'output/030_map/star-index/SA'
    output:
        ('output/030_map/{s}/'
         '{s}.Log.final.out')
    params:
        genome_dir = 'output/030_map/star-index',
        prefix = 'output/030_map/{s}/{s}.'
    threads:
        multiprocessing.cpu_count()
    log:
        'output/logs/030_map/{s}.log'
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
        'data/assemblies/flye_denovo_full.racon.fasta'
    output:
        'output/030_map/star-index/SA',
        'output/030_map/star-index/chrNameLength.txt'
    params:
        outdir = 'output/030_map/star-index'
    threads:
        multiprocessing.cpu_count()
    log:
        'output/logs/030_map/generate_genome.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.outdir} '
        '--genomeFastaFiles {input} '
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

