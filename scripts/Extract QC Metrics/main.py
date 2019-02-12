import pandas as pd
import numpy as np
import columns as cols
import utils

def main(path_read, path_write, filter_ms='QEp|QX'):

    evid = utils.read_mq_big(path_read + 'combined/txt/evidence.txt', cols.evid, '', filter_ms)
    evid = utils.fixevidence(evid)
    print("evid")

    list_raws = np.sort(np.unique(evid['raw_file']))
    list_expe = np.unique(evid['experiment'])
    n_expe = len(list_expe)
    list_expe_u, uppe_dict = utils.lowercaseexpdict(list_expe, n_expe)

    summ = utils.read_mq_big(path_read + "combined/txt/summary.txt", cols.summ, list_raws, filter_ms)
    summ_qc = utils.fixsummary(summ)
    print("summ")

    cont_qc_values, cont_qc_pct = utils.make_the_pct(evid, list_raws, 'potential_contaminant', '+', isinclude=True)
    evid = utils.dropcontaminant(evid)

    miss_qc_values, miss_qc_pct = utils.make_the_pct(evid, list_raws,  'missed_cleavages', 0)
    miss_kqc_values, miss_kqc_pct = utils.make_the_pct(evid[evid['last_aa'] == 'K'], list_raws, 'missed_cleavages', 0)
    miss_rqc_values, miss_rqc_pct = utils.make_the_pct(evid[evid['last_aa'] == 'R'], list_raws, 'missed_cleavages', 0)
    oxid_qc_values, oxid_qc_pct = utils.make_the_pct(evid, list_raws, 'oxidation_m', 0)
    type_qc_values, type_qc_pct = utils.make_the_pct(evid, list_raws, 'type', 'MULTI-MSMS', isinclude=True)

    evid_qc1 = utils.run_qc(evid, cols.evid_QC1, 'evidence.txt')
    evid_qc2 = utils.run_qc(evid[evid['type'] != 'MSMS'], cols.evid_QC2, 'evidence.txt')
    evid_qc = pd.merge(evid_qc1, evid_qc2, how='left', on='raw_file')

    mssc = utils.read_mq_big(path_read + 'combined/txt/msScans.txt', cols.mssc, list_raws, filter_ms)
    mssc = utils.fixmssc(mssc)
    mssc_qc = utils.run_qc(mssc, cols.mssc_QC, 'msScans.txt')
    print("mssc")

    msms = utils.read_mq_big(path_read + 'combined/txt/msmsScans.txt', cols.msms, list_raws, filter_ms)
    msms = utils.fixmsms(msms)


    mz_range_qc = utils.make_the_bins(msms, list_raws, 'mz_range', 'identified', '+')
    z_range_qc = utils.make_the_bins(msms, list_raws, 'charge', 'identified', '+')
    se_range_qc = utils.make_the_bins(msms, list_raws, 'scan_event_number', 'identified', '+')
    # todo: log10_range_qc = utils.make_the_bins(msms, list_raws, 'log10_tic_range', 'identified', '+')

    msms_qc = utils.run_qc(msms, cols.msms_QC, 'msmsScans.txt')
    print("msms")

    allp = utils.read_mq_big(path_read + 'combined/txt/allPeptides.txt', cols.allp, list_raws, filter_ms)
    allp = utils.fixallp(allp)
    allp_qc = utils.run_qc(allp, cols.allp_QC, "allPeptides.txt")
    print("allp")

    cols.prot = utils.protheaders(cols.prot, list_expe)
    prot = utils.read_mq_big(path_read + 'combined/txt/proteinGroups.txt', cols.prot, list_raws, filter_ms)
    prot = utils.dropcontaminant(prot)
    prot = prot.reset_index(drop=True)
    prot_qc, prot_qc2 = utils.calculateproteins(prot, list_expe_u, uppe_dict)
    print("prot")

    qc = utils.mergeandfix(summ_qc, mssc_qc, msms_qc, evid_qc, miss_qc_pct, cont_qc_pct, allp_qc, prot_qc2)
    #
    # import pdb
    # pdb.set_trace()
    qc.to_csv(path_write + 'QC_Results.tab', sep='\t', header=True, index=False)

    cont_qc_values.to_csv(path_write + 'Contaminants.tab', sep='\t', header=True, index=False)
    miss_qc_values.to_csv(path_write + 'Missed_cleavages_KR.tab', sep='\t', header=True, index=False)
    miss_kqc_values.to_csv(path_write + 'Missed_cleavages_K.tab', sep='\t', header=True, index=False)
    miss_rqc_values.to_csv(path_write + 'Missed_cleavages_R.tab', sep='\t', header=True, index=False)
    prot_qc.to_csv(path_write + 'Protein_groups.tab', sep='\t', header=True, index=False)
    type_qc_values.to_csv(path_write + 'ID_per_type.tab', sep='\t', header=True, index=False)
    mz_range_qc.to_csv(path_write + 'ID_per_mz.tab', sep='\t', header=True, index=False)
    z_range_qc.to_csv(path_write + 'ID_per_charge.tab', sep='\t', header=True, index=False)
    se_range_qc.to_csv(path_write + 'ID_per_Scan_number.tab', sep='\t', header=True, index=False)
    oxid_qc_values.to_csv(path_write + 'Oxidations.tab', sep='\t', header=True, index=False)


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description="Extract QC Metrics.")
    parser.add_argument('path_read', default='F:/Python_tests/Mq16210_ExpON_MbrON_LFQON/', help="The Folder with one MQ-analyzed RAWfile")
    parser.add_argument('path_write', default='F:/Results/', help="The Folder where we save for now the QC results in tab format")

    args = parser.parse_args()

    main(path_read=args.path_read, path_write=args.path_write)


