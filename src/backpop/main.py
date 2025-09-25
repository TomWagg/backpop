import numpy as np
from cosmic import _evolvebin
import pandas as pd

from scipy.stats import multivariate_normal
from argparse import ArgumentParser
from configparser import ConfigParser

from nautilus import Prior, Sampler

from .consts import *


class BackPop():
    def __init__(self):
        pass


def set_flags(params_in, defaults_file='cosmic_defaults.ini'):
    '''Create a Dictionary of COSMIC flags from input parameters that uses defaults from an ini file
    
    Parameters
    ----------
    params_in : dict
        Dictionary of input parameters to set. These will override the defaults
        
    Returns
    -------
    flags : dict
        Dictionary of COSMIC flags to be passed to COSMIC
    '''
    config = ConfigParser()
    config.read(defaults_file)
    flags = {section: dict(config.items(section)) for section in config.sections()}["bse"]

    natal_kick = np.zeros((2,5))
    qcrit_array = np.zeros(16)
    qc_list = ["qMSlo", "qMS", "qHG", "qGB", "qCHeB", "qAGB", "qTPAGB", "qHeMS", "qHeGB", "qHeAGB"]

    for param in params_in.keys():
        # handle kicks
        if param in ["vk1", "phi1", "theta1", "omega1", "vk2", "phi2", "theta2", "omega2"]:
            param_name = param[:-1]
            param_star = int(param[-1]) - 1
            natal_kick[param_star, NATAL_KICK_TRANSLATOR[param_name]] = params_in[param]
        elif param in qc_list:
            ind_dict = {}
            for k, v in zip(qc_list, range(0,10)):
                ind_dict[v] = k
            qcrit_array[ind_dict[param]] = params_in[param]
        else:
            flags[param] = params_in[param]


    if np.any(qcrit_array != 0.0):
        flags["qcrit_array"] = qcrit_array   
    if np.any(natal_kick != 0.0):
        flags["natal_kick_array"] = natal_kick
    
    return flags


def set_evolvebin_flags(flags):
    '''Set the flags in the _evolvebin Fortran module from a dictionary of flags
    
    Parameters
    ----------
    flags : dict
        Dictionary of COSMIC flags to be passed to COSMIC
    
    Returns
    -------
    None
    '''
    # the following is equivalent to _evolvebin.windvars.neta = flags["neta"], etc
    for g in FLAG_GROUPS:
        for k in FLAG_GROUPS[g]:
            if k not in flags:
                raise ValueError(f"flag {k} not found in flags dictionary")
            getattr(getattr(_evolvebin, g), k) = flags[k]
    return None


def select_phase(bpp, phase_select='BBH_merger'):
    '''Select the rows of the BPP array corresponding to a given phase.
    
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



def evolv2(params_in, params_out, phase_select='BBH_merger'):
    '''Evolve a binary with COSMIC given input parameters and return the output parameters
    at the time of the first BBH merger, as well as the full BPP and kick arrays.

    Parameters
    ----------
    params_in : dict
        Dictionary of input parameters that will be sampled by Nautilus
    params_out : list
        List of output parameters to return at the time of the selected phase
    phase_select : str, optional
        The phase to select the output parameters from. Currently only 'BBH_merger' is
        implemented. Default is 'BBH_merger'.
    
    Returns
    -------
    out : pd.DataFrame or None
        DataFrame of output parameters at the time of the selected phase, or None if
        the phase was not reached
    bpp : np.ndarray or None
        Full BPP array from COSMIC, or None if the phase was not reached
    kick_info : np.ndarray or None
        Full kick info array from COSMIC, or None if the phase was not reached
    '''
    # handle initial binary parameters first
    m1 = params_in["m1"] 
    m2 = params_in["m2"]
    m2, m1 = np.sort([m1,m2],axis=0)
    tb = params_in["tb"] 
    e = params_in["e"]
    # this is hardcoded for BH3... need to figure out how to specify fixed quantities..
    metallicity = 1.23e-4
    # set the other flags
    flags = set_flags(params_in)
    _ = set_evolvebin_flags(flags)
    
    bpp_columns = BPP_COLUMNS
    bcm_columns = BCM_COLUMNS
    
    col_inds_bpp = np.zeros(len(ALL_COLUMNS), dtype=int)
    col_inds_bpp[:len(bpp_columns)] = [ALL_COLUMNS.index(col) + 1 for col in bpp_columns]
    n_col_bpp = len(BPP_COLUMNS)    

    col_inds_bcm = np.zeros(len(ALL_COLUMNS), dtype=int)
    col_inds_bcm[:len(bcm_columns)] = [ALL_COLUMNS.index(col) + 1 for col in bcm_columns]
    n_col_bcm = len(BCM_COLUMNS)
    
    _evolvebin.col.n_col_bpp = n_col_bpp
    _evolvebin.col.col_inds_bpp = col_inds_bpp
    _evolvebin.col.n_col_bcm = n_col_bcm
    _evolvebin.col.col_inds_bcm = col_inds_bcm
    
    # setup the inputs for _evolvebin
    zpars = np.zeros(20)
    mass = np.array([m1,m2])
    mass0 = np.array([m1,m2])
    epoch = np.array([0.0,0.0])
    ospin = np.array([0.0,0.0])
    tphysf = 13700.0
    dtp = 0.0
    rad = np.array([0.0,0.0])
    lumin = np.array([0.0,0.0])
    massc = np.array([0.0,0.0])
    radc = np.array([0.0,0.0])
    menv = np.array([0.0,0.0])
    renv = np.array([0.0,0.0])
    B_0 = np.array([0.0,0.0])
    bacc = np.array([0.0,0.0])
    tacc = np.array([0.0,0.0])
    tms = np.array([0.0,0.0])
    bhspin = np.array([0.0,0.0])
    tphys = 0.0
    bkick = np.zeros(20)
    bpp_index_out = 0
    bcm_index_out = 0
    kick_info_out = np.zeros(34)
    kstar = np.array([1,1])
    kick_info = np.zeros((2, 18))


    [bpp_index, bcm_index, kick_info_arrays] = _evolvebin.evolv2(kstar,mass,tb,e,metallicity,tphysf,
                                                          dtp,mass0,rad,lumin,massc,radc,
                                                          menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
                                                          bhspin,tphys,zpars,bkick,kick_info)
    
    bpp = _evolvebin.binary.bpp[:bpp_index, :n_col_bpp].copy()
    _evolvebin.binary.bpp[:bpp_index, :n_col_bpp] = np.zeros(bpp.shape)
    bcm = _evolvebin.binary.bcm[:bcm_index, :n_col_bcm].copy()
    _evolvebin.binary.bcm[:bcm_index, :n_col_bcm] = np.zeros(bcm.shape)
    
    
    bpp = pd.DataFrame(bpp, columns=BPP_COLUMNS)
    
   
    kick_info = pd.DataFrame(kick_info_arrays,
                             columns=KICK_COLUMNS,
                             index=kick_info_arrays[:, -1].astype(int))
    
    out=select_phase(bpp, phase_select=phase_select)

    if len(out) > 0:
        return out[params_out].iloc[0], bpp.to_numpy(), kick_info.to_numpy()
    else:
        return None, None, None


def likelihood(rv, lower_bound, upper_bound, params_out, phase_select, x):
    '''Calculate the log-likelihood of a binary given prior bounds and input parameters
    using COSMIC to evolve the binary, select the phase of interest, and compare to
    observed binary properties

    Parameters
    ----------
    rv : scipy.stats.rv_continuous
        A scipy.stats continuous random variable object representing the likelihood
        function to evaluate the output parameters against
    lower_bound : list
        List of lower bounds for the physical parameters to enforce priors on
    upper_bound : list
        List of upper bounds for the physical parameters to enforce priors on
    params_out : list
        List of output parameters to return at the time of the selected phase
    phase_select : str, optional
        The phase to select the output parameters from. Default is 'BBH_merger'
    x : dict
        Dictionary of input parameters that will be sampled by Nautilus
    
    Returns
    -------
    ll : float
        The log-likelihood of the binary given the input parameters and priors
    bpp_flat : np.ndarray
        Flattened array of the full BPP output from COSMIC
    kick_flat : np.ndarray
        Flattened array of the full kick info output from COSMIC
    '''

    # enforce limits on physical values for kicks
    for i, name in enumerate(x):
        val = x[name]
        if name in ["theta1", "phi1", "omega1", "theta2", "phi2", "omega2"]:
            if val < lower_bound[i] or val > upper_bound[i]:
                # return invalid flattened arrays
                return -np.inf, np.full(np.prod(BPP_SHAPE), np.nan, dtype=float), np.full(np.prod(KICK_SHAPE), np.nan, dtype=float)

    # evolve the binary
    result = evolv2(x, params_out, phase_select=phase_select)
    # check result and calculate likelihood
    if result[0] is None:
        return -np.inf, np.full(np.prod(BPP_SHAPE), np.nan, dtype=float), np.full(np.prod(KICK_SHAPE), np.nan, dtype=float)
    ll = np.sum(rv.logpdf(result[0]))

    # flatten arrays and force dtype
    bpp_flat = np.array(result[1], dtype=float).ravel()
    kick_flat = np.array(result[2], dtype=float).ravel()

    # check shapes
    if bpp_flat.size != np.prod(BPP_SHAPE) or kick_flat.size != np.prod(KICK_SHAPE):
        return -np.inf, np.full(np.prod(BPP_SHAPE), np.nan, dtype=float), np.full(np.prod(KICK_SHAPE), np.nan, dtype=float)
    
    # else return the log-likelihood and flattened arrays
    return ll, bpp_flat, kick_flat


if __name__ == "__main__":
    
    # First set up the commandline args
    parser = ArgumentParser()
    parser.add_argument('--n_threads', help='pass the number of cores to use', type=int, default=1)
    parser.add_argument('--n_eff', help='pass the number of effective points for Nautilus', type=int, default=10000)
    parser.add_argument('--n_live', help='pass the number of live points for Nautilus', type=int, default=3000)
    parser.add_argument('--verbose', help='supplies verbose argument to Nautilus', type=bool, default=True)
    parser.add_argument('--filepath', help='path to output file', default="./")
    parser.add_argument('--resume', help='restart file for Nautilus', type=bool, default=False)
    args = parser.parse_args()
   

    # this should probably be an inifile...
    #Porb: 115 +/- 5
    #Eccentricity: 0.85 +/- 0.02
    #Mass 1: 9.8 +/- 2
    #Mass 2: 1 +/- 0.05
    # Initializing the multivariate normal distribution for the prior
    mean = np.array([9.8, 1, 3500.0, 0.85])
    cov = np.array([[1**2, 0, 0, 0], [0, 0.05**2, 0, 0], [0, 0, 50**2, 0], [0, 0, 0, 0.02**2]])
    rv = multivariate_normal(mean, cov)
    #m1, m2, tb, e, alpha, vk1, theta1, phi1, omega1
    m1lo = 10.0
    m2lo = 0.8
    tblo = 500.0
    elo = 0.0
    vklo = 0.0
    thetalo = 0.0
    philo = -90.0
    omegalo = 0.0
    
    m1hi = 60.0
    m2hi = 1.2
    tbhi = 5000.0
    ehi = 0.9
    vkhi = 300.0
    thetahi = 360.0
    phihi = 90.0
    omegahi = 360
    
    #m1, m2, logtb, e, alpha_1, vk1, theta1, phi1, omega1
    param_names = ["m1", "m2", "tb", "e", "vk1", "theta1", "phi1", "omega1"]
    lower_bound = np.array([m1lo, m2lo, tblo, elo, vklo, thetalo, philo, omegalo])
    upper_bound = np.array([m1hi, m2hi, tbhi, ehi, vkhi, thetahi, phihi, omegahi])
    
    resume=True
    
    phase_select="BH_MS"
    params_out = ["mass_1", "mass_2", "porb", "ecc"]
    
    # Set up Nautilus prior
    prior = Prior()
    for i in range(len(param_names)):
        prior.add_parameter(param_names[i], dist=(lower_bound[i], upper_bound[i]))
    
    num_threads = args.n_threads
    if args.verbose:
        print("using multiprocessing with " + str(num_threads) + " threads")

    # set up the blob dtype and Nautilus sampler    
    dtype = [('bpp', float, 35*len(BPP_COLUMNS)), ('kick_info', float, 2*len(KICK_COLUMNS))]
 
    sampler = Sampler(
        prior=prior, 
        likelihood=likelihood, 
        n_live=args.n_live, 
        pool=args.n_threads,
        blobs_dtype=dtype,
        filepath=args.filepath+'samples_out.hdf5', 
        resume=args.resume,
        likelihood_args=(rv, lower_bound, upper_bound, params_out, phase_select))

    sampler.run(n_eff=args.n_eff,verbose=args.verbose,discard_exploration=True)
    
    points, log_w, log_l, blobs = sampler.posterior(return_blobs=True)
    
    np.save(args.filepath+"points.npy", points)
    np.save(args.filepath+"log_w.npy", log_w)
    np.save(args.filepath+"log_l.npy", log_l)
    np.save(args.filepath+"blobs.npy", blobs)