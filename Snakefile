rule run_max_quant_analysis:
    input:
        'D:/singleRAW/{filename}/mqpar.xml',
        'D:/singleRAW/{filename}/{filename}.raw
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



