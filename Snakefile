rule copy_raw_file
    input:
        'D:/singleRAW/20190108_QX0_AlRe_MA_HeLa_500ng_LC05_1/20190108_QX0_AlRe_MA_HeLa_500ng_LC05_1.raw'
        'D:/Scripttest/20190108_QX0_AlRe_MA_HeLa_500ng_LC05_1/20190108_QX0_AlRe_MA_HeLa_500ng_LC05_1.raw'
    output:
        'D:/Scripttest/20190108_QX0_AlRe_MA_HeLa_500ng_LC05_1/20190108_QX0_AlRe_MA_HeLa_500ng_LC05_1.raw'
        'D:/Scripttest/20190108_QX0_AlRe_MA_HeLa_500ng_LC05_1'
    shell:
        "mkdir D:/Scripttest/20190108_QX0_AlRe_MA_HeLa_500ng_LC05_1"
        "copy {input[0]} {output[0]}"

rule prepare_max_quant_analysis
    input:
        "data/fasta/20190110_HomoSapiens_95965entries.fasta"
        'D:/singleRAW/{filename}/{filename}.raw'
        'D:/MaxQuant/mqpar.xml'
        analysis_directory
    params:
        threads = 2
    output:
        'D:/singleRAW/{filename}/mqpar.xml'
    shell:
        python scripts/'Prepare Max Quant/preparemaxquant.py {input[1]} {input[2]} {input3} {input4} {params.threads} {output}

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



