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
