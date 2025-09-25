__all__ = ['select_phase']

def select_phase(bpp, phase_select='BBH_merger'):
    '''Select the rows of the BPP array corresponding to a given phase.

    TODO: Allow custom masks based on different columns/operations/values
    
    Parameters
    ----------
    bpp : pd.DataFrame
        DataFrame of the BPP array from COSMIC
    phase_select : str, optional
        The phase to select. Currently only 'BBH_merger' is implemented. Default is
        'BBH_merger'.
    
    Returns
    -------
    out : pd.DataFrame or None
        DataFrame of output parameters at the time of the selected phase, or None if
        the phase was not reached
    '''
    if phase_select == 'BNS_merger':
        out = bpp.loc[(bpp.kstar_1 == 13) & (bpp.kstar_2 == 13) & (bpp.evol_type == 3)]
    elif phase_select == 'NSBH_merger':
        out = bpp.loc[((bpp.kstar_1 == 14) & (bpp.kstar_2 == 13) & (bpp.evol_type == 3)) |
                      ((bpp.kstar_1 == 13) & (bpp.kstar_2 == 14) & (bpp.evol_type == 3))]
    elif phase_select == "BBH_merger":
        out = bpp.loc[(bpp.kstar_1 == 14) & (bpp.kstar_2 == 14) & (bpp.evol_type == 3)]
    elif phase_select == "BH_MS":
        out = bpp.loc[((bpp.kstar_1 == 14) & (bpp.kstar_2.isin([0,1])) & (bpp.sep > 0)) |
                      ((bpp.kstar_1.isin([0,1])) & (bpp.kstar_2 == 14) & (bpp.sep > 0))]    
    elif phase_select == "NS_MS":
        out = bpp.loc[((bpp.kstar_1 == 13) & (bpp.kstar_2.isin([0,1])) & (bpp.sep > 0)) |
                      ((bpp.kstar_1.isin([0,1])) & (bpp.kstar_2 == 13) & (bpp.sep > 0))]
    elif phase_select == "WD_MS":
        out = bpp.loc[((bpp.kstar_1.isin([10,11,12])) & (bpp.kstar_2.isin([0,1])) & (bpp.sep > 0)) |
                      ((bpp.kstar_1.isin([0,1])) & (bpp.kstar_2.isin([10,11,12])) & (bpp.sep > 0))]
    elif phase_select == "BH_GS":
        out = bpp.loc[((bpp.kstar_1 == 14) & (bpp.kstar_2 == 3) & (bpp.sep > 0)) |
                      ((bpp.kstar_1 == 3) & (bpp.kstar_2 == 14) & (bpp.sep > 0))]
    elif phase_select == "NS_GS":
        out = bpp.loc[((bpp.kstar_1 == 13) & (bpp.kstar_2 == 3) & (bpp.sep > 0)) |
                      ((bpp.kstar_1 == 3) & (bpp.kstar_2 == 13) & (bpp.sep > 0))]
    elif phase_select == "WD_GS":
        out = bpp.loc[((bpp.kstar_1.isin([10,11,12])) & (bpp.kstar_2 == 3) & (bpp.sep > 0)) |
                      ((bpp.kstar_1 == 3) & (bpp.kstar_2.isin([10,11,12])) & (bpp.sep > 0))]
    else:
        raise ValueError("you've not specified one of the available choices for backpop stages. choose from: BNS_merger,NSBH_merger,BBH_merger,BH_MS,NS_MS,WD_MS,BH_GS,NS_GS,WD_GS")
    
    return out