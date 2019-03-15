import os
from os.path import join
import shutil
rawdirectory = 'D:\\2\\'


raw_files, =glob_wildcards(join(rawdirectory, '{filename}.raw'))
print(raw_files)

rule all:
    input: expand('D:\\2\\{filename}\\{filename}.raw', filename=raw_files),
        expand('D:\\2\\{filename}\\mqpar.xml',  filename=raw_files)


rule create_folder:
    input: expand('D:\\2\\{filename}.raw', filename=raw_files)
    output: directory(expand('D:\\2\\{filename}', filename=raw_files))
    run:
        pathdir = join(rawdirectory + raw_files[0])
        os.mkdir(pathdir)
        print (pathdir)

rule copy_raw_file:
    input:
        r'D:\\2\\'
    output:
        expand('D:\\2\\{filename}\\{filename}.raw', filename=raw_files)
    run:
        pathfrom = (join(rawdirectory + raw_files[0] + '.raw'))
        pathto = (join(rawdirectory + raw_files[0] + '\\' + raw_files[0] + '.raw'))
        shutil.copyfile(pathfrom, pathto)


rule prepare_max_quant_analysis:
    input:
        r"data/fasta/20190110_HomoSapiens_95965entries.fasta",
        r'D:/2/{filename}/{filename}.raw',
        r'C:/MQ/mqpar.xml',
        r'D:/2/{filename}'
    params:
        threads = 2
    output:
        r'D:/2/{filename}/mqpar.xml'
    shell:
        r"python scripts/prepare_max_quant_analysis/preparemaxquant.py {input[1]} {input[2]} {input[3]} {params.threads}"

rule run_max_quant_analysis:
    input:
        r'D:/1/{filename}/mqpar.xml',
        r'D:/1/{filename}/{filename}.raw',
        r"data/fasta/20190110_HomoSapiens_95965entries.fasta"
    output:
        'D:/1/{filename}/combined/txt/summary.txt'
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



