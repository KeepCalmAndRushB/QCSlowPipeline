rule extract_qc_metrics:
    input:
        'O:/20190226_QX4_JoSw_MA_HeLa_500ng_LC11/',
        'F:/Results/'
    output:
        'F:/Results/proteinGroups.tab'
    shell:
        "python scripts/'Extract QC Metrics'/main.py {input[0]} {input[1]}"




python main.py F:/Python_tests/Mq16210_ExpON_MbrON_LFQON/ F:/Results/