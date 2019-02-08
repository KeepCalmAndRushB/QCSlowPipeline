import pandas as pd
import numpy as np
import re


def evidheaders(evid_cols, isoxid, isphos, isdeam, verbose=False):
    """Depending on how the MQ analysis was done, it chooses which columns should be loaded from the evidence.txt file"""
    evid_cols = evid_cols + ['Experiment']
    evid_cols = evid_cols + ['Oxidation (M)'] if isoxid else evid_cols
    evid_cols = evid_cols + ['Phospho (STY)'] if isphos else evid_cols
    evid_cols = evid_cols + ['Deamidation (NQ)'] if isdeam else evid_cols

    if verbose:
        print('Columns that will be loaded from evidence.txt:')
        print("\n".join(evid_cols))

    return evid_cols


def protheaders(prot_cols, list_expe, isexpe, islfq, ismatc, verbose=False):
    """Depending on how the MQ analysis was done, it chooses which columns should be loaded from the proteinGroups.txt file"""
    prot_cols_isexpet = ['Intensity ' + s for s in list_expe] + ['Razor + unique peptides ' + s for s in list_expe] + ['Sequence coverage ' + s for s in list_expe + ' [%]']
    prot_cols_islfqt = ['LFQ intensity ' + s for s in list_expe]
    prot_cols_ismatct = ['Identification type ' + s for s in list_expe]
    prot_cols = prot_cols + prot_cols_isexpet if isexpe else prot_cols
    prot_cols = prot_cols + prot_cols_islfqt if islfq else prot_cols
    prot_cols = prot_cols + prot_cols_ismatct if ismatc else prot_cols

    if verbose:
        print('Columns that will be loaded from proteinGroups.txt:')
        print("\n".join(prot_cols))

    return prot_cols



def read_first(mqoutput):
    """Reads first 10 rows of a MQ txt file (all columns) changing column headers
    (some header renaming, all lowercases, no spaces, no slashes, no parenthesis)
    Used to check from evidence.txt how the MQanalysis was done (phospho? oxidations?...)"""

    df = pd.read_table(mqoutput, sep='\t', nrows=10)
    df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_').str.remove('/').str.remove('(').str.remove(')')
    return df


def read_mq_small(mqoutput, selection):
    """Reads all rows of a MQ txt file (selected columns) changing column headers
    (some header renaming, all lowercases, no spaces, no slashes, no parenthesis)"""

    df = pd.read_table(mqoutput, usecols=selection, sep='\t')
    df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_').str.remove('/').str.remove('(').str.remove(')')
    return df


def read_mq_big(mqoutput, selection, list_raws, filter_ms, time_start=0, time_end=120, boxplot_slice=10, mz_slice=100):
    """Reads all rows of a MQ txt file (selected columns) changing column headers
    (some header renaming, all lowercases, no spaces, no slashes, no parenthesis)
    If a raw_file column is present: it keeps only rows with values present in 'list_raws'
    If a raw_file column is present: it keeps only rows with instruments defined in 'filter_ms'
    If a raw_file column is present: it sorts valuesby='raw_file'
    If a retention_time column is present: it makes a 'retention_time_range' column by 'boxplot_slice' values
    If a retention_time column is present: it keeps only values >= 'time_start' and <= 'time_end'
    If present, it filters out 'reverse' and 'only_identified_by_site' 
    If a mz column is present: it makes a 'mz_range' column by 'mz_slice' values"""

    df = pd.read_table(mqoutput, usecols=selection, sep='\t')
    column_dict = {'Contaminant': 'Potential_contaminant',
                   'Av_ Absolute Mass Deviation': 'Av_ Absolute Mass Deviation [ppm]'}
    df.columns = [column_dict.get(x, x) for x in df.columns]
    df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_').str.remove('/').str.remove('(').str.remove(')')
    if ('raw_file' in list(df.columns)) & (len(list_raws) > 0):
        df = df.loc[df['raw_file'].isin(list_raws)]
        df = df.loc[df['raw_file'].str.contains(filter_ms)]
        df = df.sort_values(by='raw_file', ascending=True)
    if 'retention_time' in list(df.columns):
        df = df[(df.retention_time >= time_start) & (df.retention_time <= time_end)]
        df['retention_time_range'] = np.round(((df['retention_time'] + boxplot_slice / 2) / boxplot_slice), 0).astype(
            int) * boxplot_slice
    if 'reverse' in list(df.columns):
        df = df[df['reverse'] != '+']
        df = df.drop(columns='reverse')
    if 'only_identified_by_site' in list(df.columns):
        df = df[df['only_identified_by_site'] != '+']
        df = df.drop(columns='only_identified_by_site')
    if 'mz' in list(df.columns):
        df['mz_range'] = mz_slice * np.ceil(df['mz'] / mz_slice).astype(int)
    return df


def segm(df, time_min, time_max):
    """If a retention_time column is present: it keeps only values >= 'time_min' and <= 'time_max'"""
    if 'retention_time' in list(df.columns):
        df = df[(df.retention_time >= time_min) & (df.retention_time <= time_max)]
    return df


def make_the_bins(df, list_raws, column_values, column_of_interest, criteria, verbose=False):
    """From a dataframe (1st argument), it counts how many total values are present in the column whose header is the 4rd argument
    Then it calculates how many of these values are equal to the value given in the 5th argument.
    Calculations are split for ranges, defined in the column whose header is the 3nd argument.
    All this doesn't go to the big table, but it is great for plots"""

    qc = pd.DataFrame(columns=['rawfile', column_values, 'total', column_of_interest])
    list_values = np.sort(np.unique(df[column_values]))

    for p in range(len(list_raws)):
        for i in range(len(list_values)):
            df2 = df.copy()
            mask = (df2[column_values] == list_values[i]) & (df2['raw_file'] == list_raws[p])
            df2 = df2[mask]
            new_row = [[list_raws[p], list_values[i], len(df2), len(df2[df2[column_of_interest] == criteria])]]
            qc = qc.append(pd.DataFrame(new_row, columns=['rawfile', column_values, 'total', column_of_interest]),
                           ignore_index=True)

            if verbose:
                print(new_row)

    qc['id_rate'] = ''

    for i in range(len(qc)):
        if qc['total'][i] == 0:
            qc['id_rate'][i] = 0
        else:
            qc['id_rate'][i] = (qc[column_of_interest][i] / qc['total'][i] * 100)
    return qc


def make_the_pct(df, list_raws, column_of_interest, criteria, isinclude=False, verbose=False):

    """From a dataframe (1st argument), it counts how many total values are present in the column whose header is the 3rd argument
    Then it calculates by default how many of these values are different to the value given in the 4th argument (missed cleavages...)
    if isinclude=True, it calculates how many of these values are equal to the value given in the 4th argument (contaminants...)
    2 Outputs: [0]values, for plotting; [1]percentages: for the big table"""

    qc1 = df.groupby(['raw_file', column_of_interest]).size().to_frame('count').reset_index()
    qc2 = pd.DataFrame(columns=['raw_file', column_of_interest])

    for p in range(len(list_raws)):
        df = qc1.copy()
        df = df[df['raw_file'] == list_raws[p]]
        freq = float(df[df[column_of_interest] == criteria]['count'] / sum(df['count']) * 100)
        freq = freq if isinclude else (100 - freq)
        new_row = [[list_raws[p], freq]]
        qc2 = qc2.append(pd.DataFrame(new_row, columns=['raw_file', column_of_interest]), ignore_index=True)
        if verbose:
            print(new_row)
    return qc1, qc2


def run_qc(df, columns_of_interest, output_name, time_min=20, time_max=90):
    df = df[(df.retention_time >= time_min) & (df.retention_time <= time_max)]
    df_qc = df.groupby('raw_file')[columns_of_interest].median().reset_index()
    df_qc.columns = [
        'raw_file' if 'raw_file' in header else output_name + '_' + header + '_median_' + str(time_min) + '_' + str(
            time_max) + '_min' for header in df_qc.columns]
    return df_qc


def lowercaseexpdict(list_expe, n_expe, verbose=False):
    list_expe_u = list_expe.copy()
    for a in range(n_expe):
        list_expe_u[a] = list_expe[a].strip().lower().replace(' ', '_').replace('/', '').replace('(', '').replace(')', '')
    uppe_dict = dict(zip(list_expe_u, list_expe))

    if verbose:
        print("{" + "\n".join("{}: {}".format(k, v) for k, v in uppe_dict.items()) + "}")

    return list_expe_u, uppe_dict


def rawsexpdict(summ, verbose=False):
    expe_dict = {summ.iloc[i]['experiment']: summ.iloc[i]['raw_file'] for i in range(len(summ))}
    expe_dict = {k.strip().lower().replace(' ', '_').replace('/', '').replace('(', '').replace(')', ''): v for k, v in expe_dict.items()}
    if verbose:
        print("{" + "\n".join("{}: {}".format(k, v) for k, v in expe_dict.items()) + "}")
    return expe_dict


def fixsummary(summ):
    summ_qc = summ.copy()
    summ_qc['msms_submitted_peak/sil'] = summ_qc["msms_submitted_peak"] / summ_qc['msms_submitted_sil']
    return summ_qc


def fixevidence(evid):
    evid['resolution'] = evid['resolution'].fillna(0)
    evid['retention_length'] = evid['retention_length'] * 60
    evid['peak_Tailing_USP'] = (evid['calibrated_retention_time_finish'] - evid['calibrated_retention_time_start']) / (2 * (evid['calibrated_retention_time'] - evid['calibrated_retention_time_start']))
    evid['peak_Asymmetry_Tosoh'] = (evid['calibrated_retention_time_finish'] - evid['calibrated_retention_time']) / (evid['calibrated_retention_time'] - evid['calibrated_retention_time_start'])
    evid['potential_contaminant'] = evid['potential_contaminant'].fillna('-')
    evid['last_aa'] = evid['modified_sequence'].str[-2]
    return evid


def fixmsms(msms):
    msms['rawovftt'] = msms['rawovftt'].fillna(0)
    msms['ion_injection_time'] = msms['ion_injection_time'].fillna(0)
    msms['log10_total_ion_current'] = np.log10(msms['total_ion_current'])
    msms['log10_tic_range'] = (0.5 * np.ceil(msms['log10_total_ion_current'] / 0.5))
    return msms


def fixmssc(mssc):
    mssc['rawovftt'] = mssc['rawovftt'].fillna(0)
    mssc['ctcd_comp'] = mssc['ctcd_comp'].fillna(0)
    mssc['rawovftt_x_ctcd_comp'] = mssc['rawovftt'] * mssc['ctcd_comp']
    mssc['spray'] = (mssc['total_ion_current'] / mssc['total_ion_current'].shift(+1)).fillna(1)
    mssc['spray'] = np.where(mssc['spray'] < 1, 1 / mssc['spray'], mssc['spray'])
    mssc['spray'] = np.where(mssc['spray'] > 2, 2, mssc['spray'])
    return mssc


def fixallp(allp):
    allp['fwhm_to_base'] = allp['retention_length_fwhm'] / allp['retention_length']
    return allp


def dropcontaminant(df):
    df = df[df['potential_contaminant'] != '+']
    df = df.drop(columns='potential_contaminant')
    return df


def calculateproteins(prot, ismatc, islfq, list_expe_u, uppe_dict, verbose=False):
    if ismatc:
        type_matc = ['By MS/MS', 'By matching']
    else:
        for s in list_expe_u:
            prot['identification_type_' + s] = 'By MS/MS no MBR'
        type_matc = ['By MS/MS no MBR']

    prot_qc_columns = ['experiment', 'type', 'sequence_coverage', 'intensity', 'razor_+_unique peptides']
    prot_qc_columns = prot_qc_columns + ['lfq_intensity'] if islfq else prot_qc_columns

    if verbose:
        print('We will calculate prGroups for:')
        print(("\n".join(prot_qc_columns)) + '\n')
    
    prot_qc = pd.DataFrame(columns=prot_qc_columns)
    for p in range(len(list_expe_u)):
        for i in range(len(type_matc)):
            df = prot.copy()
            col_1 = 'sequence_coverage_' + list_expe_u[p] + '_[%]'
            col_2 = 'intensity_' + list_expe_u[p]
            col_3 = 'razor_+_unique_peptides_' + list_expe_u[p]
            col_4 = 'identification_type_' + list_expe_u[p]
            if islfq:
                col_5 = 'lfq_intensity_' + list_expe_u[p]
            sequence_coverage = np.sum((df[col_1] > 0) & (df[col_4] == type_matc[i]))
            intensity = np.sum((df[col_2] > 0) & (df[col_4] == type_matc[i]))
            razor_unique_peptides = np.sum((df[col_3] > 1) & (df[col_4] == type_matc[i]))
            if islfq:
                lfq = np.sum((df[col_5] > 0) & (df[col_4] == type_matc[i]))
            if islfq:
                new_row = [[list_expe_u[p], type_matc[i], sequence_coverage, intensity, razor_unique_peptides, lfq]]
            else:
                new_row = [[list_expe_u[p], type_matc[i], sequence_coverage, intensity, razor_unique_peptides]]
            if verbose:
                print(new_row)
            prot_qc = prot_qc.append(pd.DataFrame(new_row, columns=prot_qc_columns), ignore_index=True)
    prot_qc['experiment'].replace(uppe_dict, inplace=True)
    prot_qc2 = prot_qc.groupby(['experiment'])[
        ['sequence_coverage', 'intensity', 'razor_+_unique peptides']].sum().reset_index()
    return prot_qc, prot_qc2


def mergeandfix(summ_qc, mssc_qc, msms_qc, evid_qc, miss_qc_pct, cont_qc_pct, allp_qc, prot_qc2, file_label):
    qc = summ_qc.merge(mssc_qc, on='raw_file').merge(msms_qc, on='raw_file').merge(evid_qc, on='raw_file').merge(miss_qc_pct, on='raw_file').merge(cont_qc_pct, on='raw_file').merge(allp_qc, on='raw_file')
    qc = pd.merge(qc, prot_qc2, on='experiment')
    qc['file_label'] = file_label
    qc['comment'] = ''
    qc['instrument'] = ''
    for i in range(len(qc)):
        qc.loc[i, 'instrument'] = re.compile('QEp[0-9]|QX[0-9]|Orbi[0-9]').findall(qc['raw_file'][i])
    qc['date'] = ''
    for i in range(len(qc)):
        if bool(re.compile('\\d{12}$').search(qc['raw_file'][i])):
            qc.loc[i, 'date'] = re.compile('\\d{12}$').findall(qc['raw_file'][i])
            qc.loc[i, 'date'] = re.compile('^\\d{6}').findall(qc['date'][i])
            qc.loc[i, 'date'] = '20' + qc['date'][i]
        else:
            qc.loc[i, 'date'] = re.compile('^\\d{8}').findall(qc['raw_file'][i])

    qc = qc.iloc[:, [52, 55, 53, 54, 0, 7, 10, 2, 3, 12, 11, 4, 5, 8, 9, 6, 17, 13, 14, 16, 27, 28, 29, 32, 33, 30, 31, 34, 35, 22, 21, 18, 19, 20, 24, 23, 25, 26, 45, 44, 15, 38, 39, 40, 41, 37, 36, 42, 43, 46, 47, 48, 1, 49, 50, 51]]
    return qc


def writeresults(df, text, pre_name):
    df.to_csv(pre_name + text, sep='\t', header=True, index=False)
