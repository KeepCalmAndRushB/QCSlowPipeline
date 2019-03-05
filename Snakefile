rule copy_raw_file
    input:
        MS12
        # new file detected...
    output:
        'D:/singleRAW/{newfolder}/{newraw}'
    shell:
        "copy {input[0]} {input[1]}"

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
        python scripts/'Prepare Max Quant/preparemaxquant.py {input[1]} {input[2]} {params.threads} {output}

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
        'O:/20190228_QX1_ChDe_MA_HeLa_200ng_LC14_DMSO_2/',
        'F:/Results/'
    output:
        'F:/Results/QC_All_together.tab'
    shell:
        "python scripts/'Extract QC Metrics'/main.py {input[0]} {input[1]}"



