import os
from os.path import join, dirname
import shutil

rawdirectory = 'Y:/QHFX2/2019'
Destdir = 'Y:/2/'
maxquant_dir = 'W:/Users/Christian D/MQ16210 Maint'

raw_files, = glob_wildcards(join(rawdirectory, '{filename}.raw'))
raw_files = raw_files[:3]


rule all:
    input: expand(join(Destdir, '{filename}/added_to_database'), filename=raw_files)


rule copy_raw_file:
    input:
        join(rawdirectory, '{filename}.raw')
    output:
        join(Destdir, '{filename}/{filename}.raw')
    run:
        os.makedirs(dirname(output[0]), exist_ok=True)
        shutil.copyfile(input[0], output[0])

rule prepare_max_quant_analysis:
    input:
        join(Destdir, '{filename}/{filename}.raw'),
        join(maxquant_dir, 'mqpar.xml')
    params:
        threads = 2,
        adir = join(Destdir, '{filename}')
    output:
        join(Destdir, '{filename}/mqpar.xml'),
    shell:
        "python scripts/prepare_max_quant_analysis/preparemaxquant.py {input[0]} {input[1]} {params.threads} {params.adir}"

rule run_max_quant_analysis:
    input:
        join(Destdir, '{filename}/mqpar.xml'),
        join(Destdir, '{filename}/{filename}.raw')
    output:
        join(Destdir, '{filename}/maxquant_finished')
    shell:
        "{maxquant_dir}/bin/MaxQuantCmd.exe {input[0]} && (find /c 'Finish writing tab' D:/2/{{filename}}/combined/proc/#runningTimes.txt) && (echo yes > D:/2/{{filename}}/maxquant_finished.txt)"


rule extract_qc_metrics:
    input:
        join(Destdir, '{filename}/maxquant_finished')
    params:
        outdir = join(Destdir, '{filename}/')
    output:
        join(Destdir, '{filename}/combined/txt/QC_Results.tab')
    shell:
        "python scripts/'Extract QC Metrics'/main.py {params.outdir}"


rule add_metrics_to_database:
    input:
        join(Destdir, '{filename}/combined/txt/QC_Results.tab')
    output:
        join(Destdir, '{filename}/added_to_database')
    run:
        with open(output, 'w'):
            pass