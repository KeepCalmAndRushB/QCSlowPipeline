import os
from os.path import join
import shutil

rawdirectory = 'C:\\1\\'
Destdir = 'D:\\2\\'

raw_files, =glob_wildcards(join(rawdirectory, '{filename}.raw'))
print(raw_files)

analysis_dir = join('D:\\2\\' + raw_files[0])

rule all:
    input: expand('D:\\1\\{filename}\\combined\\txt\\summary.txt', filename=raw_files)


rule create_folder:
    input: expand('C:\\1\\{filename}.raw', filename=raw_files)
    output: directory(expand('D:\\2\\{filename}', filename=raw_files))
    run:
        pathdir = join(Destdir + raw_files[0])
        os.mkdir(pathdir)
        print (pathdir)

rule copy_raw_file:
    input:
        'D:\\2\\'
    output:
        expand('D:\\2\\{filename}\\{filename}.raw', filename=raw_files)
    run:
        pathfrom = (join(rawdirectory + raw_files[0] + '.raw'))
        pathto = (join('D:\\2\\' + raw_files[0] + '\\' + raw_files[0] + '.raw'))
        shutil.copyfile(pathfrom, pathto)

rule prepare_max_quant_analysis:
    input:
        r"data/fasta/20190110_HomoSapiens_95965entries.fasta",
        expand ('D:\\2\\{filename}\\{filename}.raw', filename=raw_files),
        r'C:/MQ/mqpar.xml',
    params:
        threads = 2,
        adir = analysis_dir
    output:
        expand('D:/2/{filename}/mqpar.xml', filename=raw_files)
    shell:
        r"python scripts/prepare_max_quant_analysis/preparemaxquant.py {input[1]} {input[2]} {params.threads} {params.adir}"

rule run_max_quant_analysis:
    input:
        expand(r'D:\2\{filename}\mqpar.xml', filename=raw_files),
        expand('D:\\2\\{filename}\\{filename}.raw', filename=raw_files),
        "data\\fasta\\20190110_HomoSapiens_95965entries.fasta"
    output:
        expand('D:\\2\\{filename}\\combined\\txt\\summary.txt', filename=raw_files)
    shell:
        "C:\\MQ\\MaxQuant\\bin\\MaxQuantCmd.exe {input[0]}"

rule extract_qc_metrics:
    input:
        'F:/Python_tests/Mq16210_ExpON_MbrON_LFQON/',
        'F:/Results/'
    output:
        'F:/Results/QC_Results.tab'
    shell:
        "python scripts/'Extract QC Metrics'/main.py {input[0]} {input[1]}"




