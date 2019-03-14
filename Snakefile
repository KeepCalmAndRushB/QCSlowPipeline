rule create_folder:
    input: r'D:\1\test.raw'
    output: directory (r'D:\1\test')
    shell:
        r'mkdir D:\1\test'

rule copy_raw_file:
    input: r'D:\1\test'
    output: r'D:\1\test\test.raw'
    shell:
        r"copy D:\1\test.raw D:\1\test\test.raw"

rule prepare_max_quant_analysis:
    input:
        r"data/fasta/20190110_HomoSapiens_95965entries.fasta",
        r'D:/1/test/test.raw',
        r'C:/MQ/mqpar.xml',
        r'D:/1/test'
    params:
        threads = 2
    output:
        r'D:/1/test/mqpar.xml'
    shell:
        r"python scripts/prepare_max_quant_analysis/preparemaxquant.py {input[1]} {input[2]} {input[3]} {params.threads}"

rule run_max_quant_analysis:
    input:
        r'D:/2/test/mqpar.xml',
        r'D:/2/test/test.raw',
        r"data/fasta/20190110_HomoSapiens_95965entries.fasta"
    output:
        'D:/Scripttest/{filename}/combined/txt/summary.txt'
    shell:
        "D:\\MaxQuant\\bin\\MaxQuantCmd.exe {input[0]}"

rule extract_qc_metrics:
    input:
        'F:/Python_tests/Mq16210_ExpON_MbrON_LFQON/',
        'F:/Results/'
    output:
        'F:/Results/QC_Results.tab'
    shell:
        "python scripts/'Extract QC Metrics'/main.py {input[0]} {input[1]}"



