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
