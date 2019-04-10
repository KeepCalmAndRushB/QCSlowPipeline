import os
from os.path import join, dirname
import shutil

configfile: "config.yaml"
workdir: 'D:/2/'


rawdir = join(config['rawdirectory'], '{filename}.raw')
raw_files, = glob_wildcards(rawdir)
raw_files = raw_files[:3]


rule all:
    input: expand('{filename}/added_to_database', filename=raw_files)


rule copy_raw_file:
    input:
        join(config['rawdirectory'], '{filename}.raw')
    output:
        '{filename}/{filename}.raw'
    run:
        os.makedirs(dirname(output), exist_ok=True)
        shutil.copyfile(input[0], output)

rule prepare_max_quant_analysis:
    input:
        '{filename}/{filename}.raw',
        join(config['maxquant_dir'], 'mqpar.xml')
    params:
        threads = 2,
        adir = '{filename}'
    output:
        '{filename}/mqpar.xml'
    shell:
        "python scripts/prepare_max_quant_analysis/preparemaxquant.py {input[0]} {input[1]} {params.threads} {params.adir}"

rule run_max_quant_analysis:
    input:
        '{filename}/mqpar.xml',
        '{filename}/{filename}.raw'
    output:
        '{filename}/maxquant_finished'
    shell:
        "{maxquant_dir}/bin/MaxQuantCmd.exe {input[0]} && (find /c 'Finish writing tab' D:/2/{{filename}}/combined/proc/#runningTimes.txt) && (echo yes > D:/2/{{filename}}/maxquant_finished.txt)"


rule extract_qc_metrics:
    input:
        '{filename}/maxquant_finished'
    params:
        outdir = '{filename}/'
    output:
        '{filename}/QC_Results.tab'
    shell:
        "python scripts/'Extract QC Metrics'/main.py {params.outdir}"


rule add_metrics_to_database:
    input:
        '{filename}/QC_Results.tab'
    output:
        '{filename}/added_to_database'
    shell:
        "python scripts/'Metrics to database'/main.py {input}"