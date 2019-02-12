rule prepare_max_quant_analysis
    input:
        "data/fasta/20190110_HomoSapiens_95965entries.fasta"
        'D:/singleRAW/{filename}/{filename}.raw'
        'D:/MaxQuant/mqpar.xml'
    params:
        threads = 2
    output:
        'D:/singleRAW/{filename}/mqpar.xml'
      shell:
       # 'modify mqpar file script'

rule run_max_quant_analysis:
    input:
        'D:/singleRAW/{filename}/mqpar.xml',
        'D:/singleRAW/{filename}/{filename}.raw',
        "data/fasta/20190110_HomoSapiens_95965entries.fasta"
    output:
        'D:/singleRAW/{filename}/combined/txt/summary.txt'
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



