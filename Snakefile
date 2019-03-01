rule extract_qc_metrics:
    input:
        'O:/20190228_QX1_ChDe_MA_HeLa_200ng_LC14_DMSO_2/',
        'F:/Results/'
    output:
        'F:/Results/QC_All_together.tab'
    shell:
        "python scripts/'Extract QC Metrics'/main.py {input[0]} {input[1]}"




python main.py F:/Python_tests/Mq16210_ExpON_MbrON_LFQON/ F:/Results/