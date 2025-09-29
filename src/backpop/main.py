import numpy as np
import pandas as pd

from scipy.stats import multivariate_normal
import os.path

from cosmic import _evolvebin
from nautilus import Prior, Sampler

from .consts import *
from .files import parse_inifile
from .phase import select_phase


class BackPop():
    """Class to sample the joint distributions of initial binary parameters and binary interaction
    assumptions using Nautilus and COSMIC.

    Parameters
    ----------
    config_file : str, optional
        Path to INI file containing configuration parameters. Default is 'params.ini'.

    Attributes
    ----------
    config_file : str
        Path to INI file containing configuration parameters.
    flags : dict
        Dictionary of COSMIC flags to be set before evolving each binary.
    config : dict
        Dictionary of backpop configuration parameters.
    obs : dict
        Dictionary of observational constraints including means, sigmas, names, and output names. Names
        correspond to the internal COSMIC fortran names, while out_names correspond to the names in the output
        BPP array.
    var : dict
        Dictionary of parameters to vary including min, max, and names. Names correspond to the COSMIC fortran
        variable names.
    fixed : dict
        Dictionary of fixed parameters and their values. Names correspond to the COSMIC fortran variable names
    rv : scipy.stats.rv_continuous
        A scipy.stats continuous random variable object representing the likelihood function to evaluate
        the output parameters against.
    prior : nautilus.Prior
        Nautilus Prior object representing the prior distributions of the parameters to be varied.
    sampler : nautilus.Sampler
        Nautilus Sampler object used to perform the sampling.
    """
    def __init__(self, config_file='params.ini'):
        self.config_file = config_file

        # parse the configuration ini file, set flags and config
        self.config, self.flags, self.obs, self.var, self.fixed = parse_inifile(self.config_file)
        self.init_flags = self.flags.copy()

        # create a scipy rv object for the likelihood
        # NOTE: currently assumes independent Gaussians (no correlated noise)
        self.rv = multivariate_normal(
            mean=np.array(self.obs["mean"]),
            cov=np.diag(np.array(self.obs["sigma"])**2)
        )
        
        # initialise the Nautilus prior
        self.prior = Prior()
        for i in range(len(self.var["name"])):
            self.prior.add_parameter(self.var["name"][i], dist=(self.var["min"][i], self.var["max"][i]))
    
    def run_sampler(self):
        """Run the Nautilus sampler to sample the joint distribution of initial binary parameters
        and binary interaction assumptions."""
        if self.config["verbose"]:
            print(f"Running sampling using multiprocessing with {self.config['n_threads']} threads")
    
        self.sampler = Sampler(
            prior=self.prior, 
            likelihood=self.likelihood, 
            n_live=self.config["n_live"], 
            pool=self.config["n_threads"],
            blobs_dtype=[('bpp', float, 35*len(BPP_COLUMNS)),
                         ('kick_info', float, 2*len(KICK_COLUMNS))],
            filepath=os.path.join(self.config["filepath"], 'samples_out.hdf5'), 
            resume=self.config["resume"]
        )
        self.sampler.run(n_eff=self.config["n_eff"], verbose=self.config["verbose"], discard_exploration=True)

        points, log_w, log_l, blobs = self.sampler.posterior(return_blobs=True)
        for label, val in zip(["points", "log_w", "log_l", "blobs"], [points, log_w, log_l, blobs]):
            np.save(os.path.join(self.config["filepath"], f"{label}.npy"), val)

    def likelihood(self, x):
        '''Calculate the log-likelihood of a binary.
        
        Calculate the log-likelihood of a binary given prior bounds and input parameters
        using COSMIC to evolve the binary, select the phase of interest, and compare to
        observed binary properties.

        Parameters
        ----------
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

        # enforce limits on physical values
        # TODO: check with Katie if this is necessary with Nautilus priors
        for i, name in enumerate(x):
            val = x[name]
            if val < self.var["min"][i] or val > self.var["max"][i]:
                # return invalid flattened arrays
                return (-np.inf, np.full(np.prod(BPP_SHAPE), np.nan, dtype=float),
                        np.full(np.prod(KICK_SHAPE), np.nan, dtype=float))

        # evolve the binary
        result = self.evolv2(x)
        # check result and calculate likelihood
        if result[0] is None:
            # print("No result!!")
            return (-np.inf, np.full(np.prod(BPP_SHAPE), np.nan, dtype=float),
                    np.full(np.prod(KICK_SHAPE), np.nan, dtype=float))
        ll = np.sum(self.rv.logpdf(result[0]))

        # flatten arrays and force dtype
        bpp_flat = np.array(result[1], dtype=float).ravel()
        kick_flat = np.array(result[2], dtype=float).ravel()

        # check shapes
        # if bpp_flat.size != np.prod(BPP_SHAPE) or kick_flat.size != np.prod(KICK_SHAPE):
        #     print(result[1].shape, result[2].shape, BPP_SHAPE, KICK_SHAPE)
        #     raise ValueError("BPP or kick array shape is incorrect")
        #     return (-np.inf, np.full(np.prod(BPP_SHAPE), np.nan, dtype=float),
        #             np.full(np.prod(KICK_SHAPE), np.nan, dtype=float))
        
        # else return the log-likelihood and flattened arrays
        return ll, bpp_flat, kick_flat
    
    def evolv2(self, params_in):
        '''Evolve a binary with COSMIC given input parameters and return the output parameters
        at the time of the first BBH merger, as well as the full BPP and kick arrays.

        Parameters
        ----------
        params_in : dict
            Dictionary of input parameters that will be sampled by Nautilus
        
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
        # handle initial binary parameters first, ensure all have been provided somewhere
        for param in ["m1", "m2", "tb", "e", "metallicity"]:
            if param not in params_in and param not in self.fixed:
                raise ValueError(f"You must provide an input value for {param} "
                                 "either as a variable or fixed parameter")
            
        # set values for the evolvebin call
        m1 = params_in["m1"] if "m1" in params_in else self.fixed["m1"]
        m2 = params_in["m2"] if "m2" in params_in else self.fixed["m2"]
        m2, m1 = np.sort([m1,m2],axis=0)
        tb = params_in["tb"] if "tb" in params_in else self.fixed["tb"]
        e = params_in["e"] if "e" in params_in else self.fixed["e"]
        metallicity = params_in["metallicity"] if "metallicity" in params_in else self.fixed["metallicity"]

        # set the other flags
        self.set_flags(params_in)
        self.set_evolvebin_flags()
        
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

        # most inputs start as a pair of zeros
        pair_vars = ["epoch", "ospin", "rad", "lumin", "massc", "radc",
                     "menv", "renv", "B_0", "bacc", "tacc", "tms", "bhspin"]
        p = {var: np.zeros(2) for var in pair_vars}

        # masses and kstars are actually set
        p["mass"] = np.array([m1, m2])
        p["mass0"] = np.array([m1, m2])
        p["kstar"] = np.array([1, 1])
        
        # setup the inputs for _evolvebin
        zpars = np.zeros(20)
        tphysf = 13700.0
        dtp = 0.0
        tphys = 0.0
        bkick = np.zeros(20)
        kick_info = np.zeros((2, 18))

        # run COSMIC!
        [bpp_index, bcm_index, kick_info_arrays] = _evolvebin.evolv2(p["kstar"], p["mass"], tb, e,
                                                                     metallicity, tphysf, dtp, p["mass0"],
                                                                     p["rad"], p["lumin"], p["massc"],
                                                                     p["radc"], p["menv"], p["renv"],
                                                                     p["ospin"], p["B_0"], p["bacc"],
                                                                     p["tacc"], p["epoch"], p["tms"],
                                                                     p["bhspin"], tphys, zpars,
                                                                     bkick, kick_info)

        bpp = _evolvebin.binary.bpp[:35, :n_col_bpp].copy()
        _evolvebin.binary.bpp[:bpp_index, :n_col_bpp] = np.zeros((bpp_index, n_col_bpp))
        bcm = _evolvebin.binary.bcm[:35, :n_col_bcm].copy()
        _evolvebin.binary.bcm[:bcm_index, :n_col_bcm] = np.zeros((bcm_index, n_col_bcm))
        # print(bpp.shape)

        bpp = pd.DataFrame(bpp, columns=BPP_COLUMNS)

        kick_info = pd.DataFrame(kick_info_arrays,
                                columns=KICK_COLUMNS,
                                index=kick_info_arrays[:, -1].astype(int))

        out = select_phase(bpp, condition=self.config["phase_condition"])

        if len(out) > 0:
            # print("FOUND ONE!")
            return out[self.obs["out_name"]].iloc[0], bpp.to_numpy(), kick_info.to_numpy()
        else:
            return None, None, None


    def set_flags(self, params_in):
        '''update dictionary of COSMIC flags with input parameters
        '''
        natal_kick = np.zeros((2,5))
        qcrit_array = np.zeros(16)
        qc_list = ["qMSlo", "qMS", "qHG", "qGB", "qCHeB", "qAGB", "qTPAGB", "qHeMS", "qHeGB", "qHeAGB"]

        # update flags based on input params
        for param in params_in.keys():
            # create natal kick arrays for each star if necessary
            if param in ["vk1", "phi1", "theta1", "omega1", "vk2", "phi2", "theta2", "omega2"]:
                param_name = param[:-1]
                param_star = int(param[-1]) - 1
                natal_kick[param_star, NATAL_KICK_TRANSLATOR[param_name]] = params_in[param]
            # same for qcrit_arrays
            elif param in qc_list:
                ind_dict = {}
                for k, v in zip(qc_list, range(0,10)):
                    ind_dict[v] = k
                qcrit_array[ind_dict[param]] = params_in[param]
            # otherwise just set the flag
            else:
                self.flags[param] = params_in[param]

        # if we set any of the arrays, update the flags
        if np.any(qcrit_array != 0.0):
            self.flags["qcrit_array"] = qcrit_array   
        if np.any(natal_kick != 0.0):
            self.flags["natal_kick_array"] = natal_kick


    def set_evolvebin_flags(self):
        '''Set the flags in the _evolvebin Fortran module from a dictionary of flags
        
        Parameters
        ----------
        flags : dict
            Dictionary of COSMIC flags to be passed to COSMIC
        '''
        # the following is equivalent to _evolvebin.windvars.neta = flags["neta"], etc
        for g in FLAG_GROUPS:
            for k in FLAG_GROUPS[g]:
                if k not in self.flags:
                    raise ValueError(f"flag {k} not found in flags dictionary")
                setattr(getattr(_evolvebin, g), k, self.flags[k])
        return None
