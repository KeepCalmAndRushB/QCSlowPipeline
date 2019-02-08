import pandas as pd
import numpy as np
import utils
import headers


time_start = 0
time_min = 20
time_max = 90
time_end = 120
boxplot_slice = 10
mz_slice = 100

path_read = 'F:/Python_tests/Mq16210_ExpON_MbrON_LFQON/'
# path_read = 'O:/20190118_QX4_AnBr_MA_HeLa_500ng_LC06_190119103719/'

path_write = 'F:/Results/'
file_label = 'ThatthewayIlikeit'
filter_MS = 'QEp|QX|Orbi'

runn_cols = ['Job', 'Running time [min]']
para_cols = ['Parameter', 'Value']
summ_cols = ['Raw file', 'Experiment', 'MS', 'MS/MS', 'MS/MS Submitted (SIL)', 'MS/MS Submitted (PEAK)',
             'MS/MS Identified [%]', 'MS/MS Identified (SIL) [%]', 'MS/MS Identified (PEAK) [%]',
             'MS/MS Identified', 'Peptide Sequences Identified', 'Peaks', 'Isotope Patterns',
             'Isotope Patterns Sequenced (z>1)', 'Isotope Patterns Sequenced (z>1) [%]',
             'Isotope Patterns Repeatedly Sequenced [%]',
             'Av. Absolute Mass Deviation [ppm]']
mssc_cols = ['Raw file', 'Retention time', 'Ion injection time', 'Total ion current', 'MS/MS count', 'Cycle time',
             'AGC Fill', 'Base peak intensity', 'RawOvFtT', 'Isotope patterns / s', 'MS/MS identification rate [%]',
             'CTCD Comp']
msms_cols = ['Raw file', 'Retention time', 'Ion injection time', 'Total ion current', 'Filtered peaks',
             'Parent intensity fraction', 'Fraction of total spectrum', 'AGC Fill', 'RawOvFtT', 'm/z', 'Charge',
             'Identified', 'Base peak intensity', 'Scan event number', 'Precursor apex fraction',
             'Precursor apex offset', 'Precursor apex offset time', 'Precursor intensity']
evid_cols = ['Raw file', 'Retention time', 'Reverse', 'Potential contaminant', 'Missed cleavages', 'Retention length',
             'Resolution', 'Mass error [ppm]',
             'Score', 'Delta score', 'Intensity', 'Type', 'Uncalibrated mass error [ppm]', 'm/z',
             'Calibrated retention time', 'Calibrated retention time start', 'Calibrated retention time finish',
             'Modified sequence', 'Charge',
             'Retention time calibration', 'MS/MS count']
allp_cols = ['Raw file', 'Retention time', 'Number of scans', 'Retention length', 'Type', 'Retention length (FWHM)',
             'Intensity', 'MSMS Scan Numbers', 'Sequence', 'Charge', 'm/z']
prot_cols = ['Protein IDs', 'Protein names', 'Gene names', 'Fasta headers', 'Number of proteins', 'Mol. weight [kDa]',
             'Sequence lengths', 'Q-value', 'Score', 'MS/MS count', 'Only identified by site',
             'Reverse', 'Potential contaminant']

file_dict = {
    'para': 'parameters.txt',
    'mzra': 'mzRange.txt',
    'summ': 'summary.txt',
    'mssc': 'msScans.txt',
    'msms': 'msmsScans.txt',
    'evid': 'evidence.txt',
    'allp': 'allPeptides.txt',
    'prot': 'proteinGroups.txt'
}

# print("{" + "\n".join("{}: {}".format(k, v) for k, v in file_dict.items()) + "}")
# print("\n")

evid_QCcols1 = ['retention_length', 'resolution', 'uncalibrated_mass_error_[ppm]', 'score', 'delta_score', 'intensity']
evid_QCcols2 = ['peak_Tailing_USP', 'peak_Asymmetry_Tosoh']
mssc_QCcols = ['ion_injection_time', 'total_ion_current', 'base_peak_intensity', 'msms_count', 'cycle_time', 'agc_fill', 'rawovftt_x_ctcd_comp', 'ctcd_comp', 'spray']
msms_QCcols = ['ion_injection_time', 'total_ion_current', 'base_peak_intensity', 'filtered_peaks', 'parent_intensity_fraction', 'rawovftt', 'agc_fill', 'precursor_apex_fraction', 'fraction_of_total_spectrum']
allp_QCcols = ['number_of_scans', 'retention_length', 'fwhm_to_base']


runn = utils.read_mq_small(path_read + 'combined/proc/#runningTimes.txt', runn_cols)

para = utils.read_mq_small(path_read + 'combined/txt/parameters.txt', para_cols)

evid = utils.read_first(path_read + 'combined/txt/evidence.txt')

isMatc = list(para['value'][para['parameter'] == 'Match between runs'] == 'True')[0]
mq_vers = list(para['value'][para['parameter'] == 'Version'])[0]

isOxid = True if 'oxidation_m' in list(evid.columns) else False
isPhos = True if 'phospho_STY' in list(evid.columns) else False
isDeam = True if 'deamidation_NQ' in list(evid.columns) else False

evid_cols = utils.evidheaders(evid_cols, isOxid, isPhos, isDeam)

print('loading evidence.txt\n')
evid = utils.read_mq_big(path_read + 'combined/txt/evidence.txt', evid_cols, '', filter_MS, time_start, time_end, boxplot_slice, mz_slice)
evid = utils.fixevidence(evid)

list_raws = np.sort(np.unique(evid['raw_file']))
n_raws = len(list_raws)

list_expe = np.unique(evid['experiment'])
n_expe = len(list_expe)

list_type = np.unique(evid['type'])
list_char = np.unique(evid['charge'])

list_expe_u, uppe_dict = utils.lowercaseexpdict(list_expe, n_expe)

print('loading summary.txt')
summ = utils.read_mq_big(path_read + "combined/txt/summary.txt", summ_cols, list_raws, filter_MS, time_start, time_end, boxplot_slice, mz_slice)
summ_qc = utils.fixsummary(summ)
print('summary.txt done\n')

expe_dict = utils.rawsexpdict(summ)

print('Little summary:')
print('Files will be read from here ............ {}'.format(path_read))
print('Results will be written here from here .. {}'.format(path_write))
print('Median calculations here ................ {} _____ {} ============ {} _____ {} min'.format(time_start, time_min, time_max, time_end))
print('We consider only the following MS ....... {}'.format(filter_MS))
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

evid_qc1 = utils.run_qc(evid, evid_QCcols1, 'evidence.txt', time_min, time_max)
evid_qc2 = utils.run_qc(evid[evid['type'] != 'MSMS'], evid_QCcols2, 'evidence.txt', time_min, time_max)
evid_qc = pd.merge(evid_qc1, evid_qc2, how='left', on='raw_file')
print('evidence.txt done\n')

print('loading msScans.txt')
mssc = utils.read_mq_big(path_read + 'combined/txt/msScans.txt', mssc_cols, list_raws, filter_MS, time_start, time_end, boxplot_slice, mz_slice)
mssc = utils.fixmssc(mssc)
mssc_qc = utils.run_qc(mssc, mssc_QCcols, 'msScans.txt', time_min, time_max)
print('msScans.txt done\n')

print('loading msmsScans.txt')
msms = utils.read_mq_big(path_read + 'combined/txt/msmsScans.txt', msms_cols, list_raws, filter_MS, time_start, time_end, boxplot_slice, mz_slice)
msms = utils.fixmsms(msms)

print('Calulating % ID per mz range')
mz_range_qc = utils.make_the_bins(msms, list_raws, 'mz_range', 'identified', '+')

print('Calulating % ID per charge state')
z_range_qc = utils.make_the_bins(msms, list_raws, 'charge', 'identified', '+')

print('Calulating % ID per scan event')
se_range_qc = utils.make_the_bins(msms, list_raws, 'scan_event_number', 'identified', '+')

print('Calulating % ID per log10 Intensity')   # todo: this calculation
# log10_range_qc = utils.make_the_bins(msms, list_raws, 'log10_tic_range', 'identified', '+')

msms_qc = utils.run_qc(msms, msms_QCcols, 'msmsScans.txt', time_min, time_max)
print('msmsScans.txt done\n')

print('loading allPeptides.txt')
allp = utils.read_mq_big(path_read + 'combined/txt/allPeptides.txt', allp_cols, list_raws, filter_MS, time_start, time_end, boxplot_slice, mz_slice)
allp = utils.fixallp(allp)
allp_qc = utils.run_qc(allp, allp_QCcols, "allPeptides.txt", time_min, time_max)
print('allPeptides.txt done\n')

print('loading proteinGroups.txt')
prot_cols = utils.protheaders(prot_cols, list_expe, isExpe, isLfq, isMatc)
prot = utils.read_mq_big(path_read + 'combined/txt/proteinGroups.txt', prot_cols, list_raws, filter_MS, time_start, time_end, boxplot_slice, mz_slice)
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
