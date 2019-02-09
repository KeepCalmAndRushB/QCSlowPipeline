runn = ['Job', 'Running time [min]']
para = ['Parameter', 'Value']
summ = ['Raw file', 'Experiment', 'MS', 'MS/MS', 'MS/MS Submitted (SIL)', 'MS/MS Submitted (PEAK)',
             'MS/MS Identified [%]', 'MS/MS Identified (SIL) [%]', 'MS/MS Identified (PEAK) [%]',
             'MS/MS Identified', 'Peptide Sequences Identified', 'Peaks', 'Isotope Patterns',
             'Isotope Patterns Sequenced (z>1)', 'Isotope Patterns Sequenced (z>1) [%]',
             'Isotope Patterns Repeatedly Sequenced [%]',
             'Av. Absolute Mass Deviation [ppm]']
mssc = ['Raw file', 'Retention time', 'Ion injection time', 'Total ion current', 'MS/MS count', 'Cycle time',
             'AGC Fill', 'Base peak intensity', 'RawOvFtT', 'Isotope patterns / s', 'MS/MS identification rate [%]',
             'CTCD Comp']
msms = ['Raw file', 'Retention time', 'Ion injection time', 'Total ion current', 'Filtered peaks',
             'Parent intensity fraction', 'Fraction of total spectrum', 'AGC Fill', 'RawOvFtT', 'm/z', 'Charge',
             'Identified', 'Base peak intensity', 'Scan event number', 'Precursor apex fraction',
             'Precursor apex offset', 'Precursor apex offset time', 'Precursor intensity']
evid = ['Raw file', 'Retention time', 'Reverse', 'Potential contaminant', 'Missed cleavages', 'Retention length',
             'Resolution', 'Mass error [ppm]',
             'Score', 'Delta score', 'Intensity', 'Type', 'Uncalibrated mass error [ppm]', 'm/z',
             'Calibrated retention time', 'Calibrated retention time start', 'Calibrated retention time finish',
             'Modified sequence', 'Charge',
             'Retention time calibration', 'MS/MS count', 'Experiment', 'Oxidation (M)']
allp = ['Raw file', 'Retention time', 'Number of scans', 'Retention length', 'Type', 'Retention length (FWHM)',
             'Intensity', 'MSMS Scan Numbers', 'Sequence', 'Charge', 'm/z']
prot = ['Protein IDs', 'Protein names', 'Gene names', 'Fasta headers', 'Number of proteins', 'Mol. weight [kDa]',
             'Sequence lengths', 'Q-value', 'Score', 'MS/MS count', 'Only identified by site',
             'Reverse', 'Potential contaminant']

evid_QC1 = ['retention_length', 'resolution', 'uncalibrated_mass_error_[ppm]', 'score', 'delta_score', 'intensity']
evid_QC2 = ['peak_Tailing_USP', 'peak_Asymmetry_Tosoh']
mssc_QC = ['ion_injection_time', 'total_ion_current', 'base_peak_intensity', 'msms_count', 'cycle_time', 'agc_fill', 'rawovftt_x_ctcd_comp', 'ctcd_comp', 'spray']
msms_QC = ['ion_injection_time', 'total_ion_current', 'base_peak_intensity', 'filtered_peaks', 'parent_intensity_fraction', 'rawovftt', 'agc_fill', 'precursor_apex_fraction', 'fraction_of_total_spectrum']
allp_QC = ['number_of_scans', 'retention_length', 'fwhm_to_base']