import os
from os.path import join, dirname
import shutil

rawdirectory = 'C:/1/'
Destdir = 'D:/2/'

raw_files = glob_wildcards(join(rawdirectory, '{filename}.raw'))
print(raw_files)


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
        'C:/MQ/mqpar.xml'
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
        "C:/MQ/MaxQuant/bin/MaxQuantCmd.exe {input[0]}"
        "find /c "Finish writing tab" D:/2/{filename}/combined/proc/#runningTimes.txt && (echo yes > D:/2/{filename}/maxquant_finished.txt)"

rule extract_qc_metrics:
    input:
        join(Destdir, '{filename}/'),
        join(Destdir, '{filename}/maxquant_finished')
    output:
        join(Destdir, '{filename}/combined/txt/QC_Results.tab')
    shell:
        "python scripts/'Extract QC Metrics'/main.py {input[0]}"


rule add_metrics_to_database:
    input:
        join(Destdir, '{filename}/combined/txt/QC_Results.tab')
    output:
        join(Destdir, '{filename}/added_to_database')
    run:
        with open(output, 'w'):
            pass

rule all:
    input: expand(join(Destdir, '{filename}/added_to_database'), filename=raw_files)