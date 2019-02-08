import pandas as pd
import numpy as np
import columns as cols
import utils

path_read = 'F:/Python_tests/Mq16210_ExpON_MbrON_LFQON/'
path_write = 'F:/Results/'
filter_MS = 'QEp|QX'




runn = utils.read_mq_small(path_read + 'combined/proc/#runningTimes.txt', cols.runn)

para = utils.read_mq_small(path_read + 'combined/txt/parameters.txt', cols.para)

evid = utils.read_first(path_read + 'combined/txt/evidence.txt')

isMatc = list(para['value'][para['parameter'] == 'Match between runs'] == 'True')[0]
mq_vers = list(para['value'][para['parameter'] == 'Version'])[0]

isOxid = True if 'oxidation_m' in list(evid.columns) else False
isPhos = True if 'phospho_STY' in list(evid.columns) else False
isDeam = True if 'deamidation_NQ' in list(evid.columns) else False

evid = utils.evidheaders(cols.evid, isOxid, isPhos, isDeam)

print('loading evidence.txt\n')
evid = utils.read_mq_big(path_read + 'combined/txt/evidence.txt', cols.evid, '', filter_MS)
evid = utils.fixevidence(evid)

list_raws = np.sort(np.unique(evid['raw_file']))
n_raws = len(list_raws)

list_expe = np.unique(evid['experiment'])
n_expe = len(list_expe)

list_type = np.unique(evid['type'])
list_char = np.unique(evid['charge'])

list_expe_u, uppe_dict = utils.lowercaseexpdict(list_expe, n_expe)

print('loading summary.txt')
summ = utils.read_mq_big(path_read + "combined/txt/summary.txt", cols.summ, list_raws, filter_MS)
summ_qc = utils.fixsummary(summ)
print('summary.txt done\n')

expe_dict = utils.rawsexpdict(summ)

print('Calulating % contaminant')
cont_qc_values, cont_qc_pct = utils.make_the_pct(evid, list_raws, 'potential_contaminant', '+', isinclude=True)
evid = utils.dropcontaminant(evid)

print('Calulating % missed_cleavages')
miss_qc_values, miss_qc_pct = utils.make_the_pct(evid, list_raws,  'missed_cleavages', 0)

print('Calulating % missed_cleavages for K')
miss_kqc_values, miss_kqc_pct = utils.make_the_pct(evid[evid['last_aa'] == 'K'], list_raws, 'missed_cleavages', 0)

print('Calulating % missed_cleavages for R')
miss_rqc_values, miss_rqc_pct = utils.make_the_pct(evid[evid['last_aa'] == 'R'], list_raws, 'missed_cleavages', 0)

if isPhos:
    print('Calulating % Phospho enrichment')
    phos_qc_values, phos_qc_pct = utils.make_the_pct(evid, list_raws, 'phospho_sty', 0)

if isDeam:
    print('Calulating % deamidation_nq')
    deam_qc_values, deam_qc_pct = utils.make_the_pct(evid, list_raws, 'deamidation_nq', 0)

if isOxid:
    print('Calulating % oxidation_m')
    oxid_qc_values, oxid_qc_pct = utils.make_the_pct(evid, list_raws, 'oxidation_m', 0)

print('Calulating % ID type\n')
type_qc_values, type_qc_pct = utils.make_the_pct(evid, list_raws, 'type', 'MULTI-MSMS', isinclude=True)

evid_qc1 = utils.run_qc(evid, cols.evid1, 'evidence.txt')
evid_qc2 = utils.run_qc(evid[evid['type'] != 'MSMS'], cols.evid2, 'evidence.txt')
evid_qc = pd.merge(evid_qc1, evid_qc2, how='left', on='raw_file')
print('evidence.txt done\n')

print('loading msScans.txt')
mssc = utils.read_mq_big(path_read + 'combined/txt/msScans.txt', cols.mssc, list_raws, filter_MS)
mssc = utils.fixmssc(mssc)
mssc_qc = utils.run_qc(mssc, cols.mssc, 'msScans.txt')
print('msScans.txt done\n')

print('loading msmsScans.txt')
msms = utils.read_mq_big(path_read + 'combined/txt/msmsScans.txt', cols.msms, list_raws, filter_MS)
msms = utils.fixmsms(msms)

print('Calulating % ID per mz range')
mz_range_qc = utils.make_the_bins(msms, list_raws, 'mz_range', 'identified', '+')

print('Calulating % ID per charge state')
z_range_qc = utils.make_the_bins(msms, list_raws, 'charge', 'identified', '+')

print('Calulating % ID per scan event')
se_range_qc = utils.make_the_bins(msms, list_raws, 'scan_event_number', 'identified', '+')

print('Calulating % ID per log10 Intensity')   # todo: this calculation
# log10_range_qc = utils.make_the_bins(msms, list_raws, 'log10_tic_range', 'identified', '+')

msms_qc = utils.run_qc(msms, cols.msms, 'msmsScans.txt')
print('msmsScans.txt done\n')

print('loading allPeptides.txt')
allp = utils.read_mq_big(path_read + 'combined/txt/allPeptides.txt', cols.allp, list_raws, filter_MS)
allp = utils.fixallp(allp)
allp_qc = utils.run_qc(allp, cols.allp, "allPeptides.txt")
print('allPeptides.txt done\n')

print('loading proteinGroups.txt')
prot = utils.protheaders(prot, list_expe, isExpe, isLfq, isMatc)
prot = utils.read_mq_big(path_read + 'combined/txt/proteinGroups.txt', cols.prot, list_raws, filter_MS)
prot = utils.dropcontaminant(prot)
prot = prot.reset_index(drop=True)
prot_qc, prot_qc2 = utils.calculateproteins(prot, isMatc, isLfq, list_expe_u, uppe_dict)
print('proteinGroups.txt done')

# todo: fix export
qc = utils.mergeandfix(summ_qc, mssc_qc, msms_qc, evid_qc, miss_qc_pct, cont_qc_pct, allp_qc, prot_qc2, file_label)



# utils.writeresults(qc, '_QC_Results.tab', pre_name)
# utils.writeresults(cont_qc_values, 'Contaminants.tab', pre_name)
# utils.writeresults(miss_qc_values, 'Missed_cleavages_KR.tab', pre_name)
# utils.writeresults(miss_kqc_values, 'Missed_cleavages_K.tab', pre_name)
# utils.writeresults(miss_rqc_values, 'Missed_cleavages_R.tab', pre_name)
# utils.writeresults(prot_qc, 'Protein_groups.tab', pre_name)
# utils.writeresults(type_qc_values, 'ID_per_type.tab', pre_name)
# utils.writeresults(mz_range_qc, 'ID_per_mz.tab', pre_name)
# utils.writeresults(z_range_qc, 'ID_per_charge.tab', pre_name)
# utils.writeresults(se_range_qc, 'ID_per_Scan_number.tab', pre_name)
# if isPhos:
#     utils.writeresults(phos_qc_values, 'Phospho.tab', pre_name)
# if isDeam:
#     utils.writeresults(deam_qc_values, 'Deamidations.tab', pre_name)
# if isOxid:
#     utils.writeresults(oxid_qc_values, 'Oxidations.tab', pre_name)
#
