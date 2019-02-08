
rule extract_qc_metrics:
    input:
        'F:/Python_tests/Mq16210_ExpON_MbrON_LFQON/',
        'F:/Results/'
    output:
        'F:/Results/QC_Results.tab'
    shell:
        "python scripts/'Extract QC Metrics'/main.py {input[0]} {input[1]}"