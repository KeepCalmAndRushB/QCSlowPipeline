import time
import pandas as pd
import numpy as np
import utils
import checks
import headers
import params as p
time_01 = time.time()

checks.directory(p.path_read)
checks.directory(p.path_write)

runn = utils.read_mq_small(p.path_read + 'combined/proc/#runningTimes.txt', p.runn_cols)
isRunn = checks.mqdone(runn, p.path_read, start_delay=2)

para = utils.read_mq_small(p.path_read + 'combined/txt/parameters.txt', p.para_cols)
isPeak = checks.peakproperties(para)

evid = utils.read_first(p.path_read + 'combined/txt/evidence.txt')
isExpe = checks.experiments(evid)

isLfq = checks.lfqison(p.path_read + 'mqpar.xml')

isMatc = list(para['value'][para['parameter'] == 'Match between runs'] == 'True')[0]
mq_vers = list(para['value'][para['parameter'] == 'Version'])[0]

isOxid = True if 'oxidation_m' in list(evid.columns) else False
isPhos = True if 'phospho_STY' in list(evid.columns) else False
isDeam = True if 'deamidation_NQ' in list(evid.columns) else False

p.evid_cols = headers.evidheaders(p.evid_cols, isOxid, isPhos, isDeam)

print('loading evidence.txt\n')
evid = utils.read_mq_big(p.path_read + 'combined/txt/evidence.txt', p.evid_cols, '', p.filter_MS, p.time_start, p.time_end, p.boxplot_slice, p.mz_slice)
evid = utils.fixevidence(evid)

list_raws = np.sort(np.unique(evid['raw_file']))
n_raws = len(list_raws)

list_expe = np.unique(evid['experiment'])
n_expe = len(list_expe)

list_type = np.unique(evid['type'])
list_char = np.unique(evid['charge'])

checks.nrawsandexper(list_raws, list_expe)

list_expe_u, uppe_dict = utils.lowercaseexpdict(list_expe, n_expe)

print('loading summary.txt')
summ = utils.read_mq_big(p.path_read + "combined/txt/summary.txt", p.summ_cols, list_raws, p.filter_MS, p.time_start, p.time_end, p.boxplot_slice, p.mz_slice)
summ_qc = utils.fixsummary(summ)
print('summary.txt done\n')

expe_dict = utils.rawsexpdict(summ)
checks.uniqueexper(expe_dict, n_raws)

print('Little summary:')
print('Files will be read from here ............ {}'.format(p.path_read))
print('Results will be written here from here .. {}'.format(p.path_write))
print('Median calculations here ................ {} _____ {} ============ {} _____ {} min'.format(p.time_start, p.time_min, p.time_max, p.time_end))
print('We consider only the following MS ....... {}'.format(p.filter_MS))
print('MQ version .............................. {}'.format(mq_vers))
print('Experiment* ............................. {}'.format(isExpe))
print('Calculate peak properties* .............. {}'.format(isPeak))
print('Match between runs ...................... {}'.format(isMatc))
print('LFQ is on ............................... {}'.format(isLfq))
print('Oxidation (M) ........................... {}'.format(isOxid))
print('Phospho (STY) ........................... {}'.format(isPhos))
print('Deamidation (NQ) ........................ {}\n'.format(isDeam))

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

evid_qc1 = utils.run_qc(evid, p.evid_QCcols1, 'evidence.txt', p.time_min, p.time_max)
evid_qc2 = utils.run_qc(evid[evid['type'] != 'MSMS'], p.evid_QCcols2, 'evidence.txt', p.time_min, p.time_max)
evid_qc = pd.merge(evid_qc1, evid_qc2, how='left', on='raw_file')
print('evidence.txt done\n')

print('loading msScans.txt')
mssc = utils.read_mq_big(p.path_read + 'combined/txt/msScans.txt', p.mssc_cols, list_raws, p.filter_MS, p.time_start, p.time_end, p.boxplot_slice, p.mz_slice)
mssc = utils.fixmssc(mssc)
mssc_qc = utils.run_qc(mssc, p.mssc_QCcols, 'msScans.txt', p.time_min, p.time_max)
print('msScans.txt done\n')

print('loading msmsScans.txt')
msms = utils.read_mq_big(p.path_read + 'combined/txt/msmsScans.txt', p.msms_cols, list_raws, p.filter_MS, p.time_start, p.time_end, p.boxplot_slice, p.mz_slice)
msms = utils.fixmsms(msms)

print('Calulating % ID per mz range')
mz_range_qc = utils.make_the_bins(msms, list_raws, 'mz_range', 'identified', '+')

print('Calulating % ID per charge state')
z_range_qc = utils.make_the_bins(msms, list_raws, 'charge', 'identified', '+')

print('Calulating % ID per scan event')
se_range_qc = utils.make_the_bins(msms, list_raws, 'scan_event_number', 'identified', '+')

print('Calulating % ID per log10 Intensity')   # todo: this calculation
# log10_range_qc = utils.make_the_bins(msms, list_raws, 'log10_tic_range', 'identified', '+')

msms_qc = utils.run_qc(msms, p.msms_QCcols, 'msmsScans.txt', p.time_min, p.time_max)
print('msmsScans.txt done\n')

print('loading allPeptides.txt')
allp = utils.read_mq_big(p.path_read + 'combined/txt/allPeptides.txt', p.allp_cols, list_raws, p.filter_MS, p.time_start, p.time_end, p.boxplot_slice, p.mz_slice)
allp = utils.fixallp(allp)
allp_qc = utils.run_qc(allp, p.allp_QCcols, "allPeptides.txt", p.time_min, p.time_max)
print('allPeptides.txt done\n')

print('loading proteinGroups.txt')
p.prot_cols = headers.protheaders(p.prot_cols, list_expe, isExpe, isLfq, isMatc)
prot = utils.read_mq_big(p.path_read + 'combined/txt/proteinGroups.txt', p.prot_cols, list_raws, p.filter_MS, p.time_start, p.time_end, p.boxplot_slice, p.mz_slice)
prot = utils.dropcontaminant(prot)
prot = prot.reset_index(drop=True)
prot_qc, prot_qc2 = utils.calculateproteins(prot, isMatc, isLfq, list_expe_u, uppe_dict)
print('proteinGroups.txt done')

qc = utils.mergeandfix(summ_qc, mssc_qc, msms_qc, evid_qc, miss_qc_pct, cont_qc_pct, allp_qc, prot_qc2, p.file_label)
pre_name = checks.writingdirectory(p.path_write, p.file_label)

utils.writeresults(qc, '_QC_Results.tab', pre_name)
utils.writeresults(cont_qc_values, 'Contaminants.tab', pre_name)
utils.writeresults(miss_qc_values, 'Missed_cleavages_KR.tab', pre_name)
utils.writeresults(miss_kqc_values, 'Missed_cleavages_K.tab', pre_name)
utils.writeresults(miss_rqc_values, 'Missed_cleavages_R.tab', pre_name)
utils.writeresults(prot_qc, 'Protein_groups.tab', pre_name)
utils.writeresults(type_qc_values, 'ID_per_type.tab', pre_name)
utils.writeresults(mz_range_qc, 'ID_per_mz.tab', pre_name)
utils.writeresults(z_range_qc, 'ID_per_charge.tab', pre_name)
utils.writeresults(se_range_qc, 'ID_per_Scan_number.tab', pre_name)
if isPhos:
    utils.writeresults(phos_qc_values, 'Phospho.tab', pre_name)
if isDeam:
    utils.writeresults(deam_qc_values, 'Deamidations.tab', pre_name)
if isOxid:
    utils.writeresults(oxid_qc_values, 'Oxidations.tab', pre_name)

print('\nAll done in {} seconds'.format(round((time.time() - time_01), 1)))
print('You can find all the tables here: ' + p.path_write + p.file_label + '/')

# todo: in the MQ txt folder create a file "Python_QC_done.tab"