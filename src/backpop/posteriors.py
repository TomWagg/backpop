import corner
import matplotlib.pyplot as plt
import numpy as np
import h5py as h5
import pandas as pd
from .consts import BPP_SHAPE, KICK_SHAPE, BPP_COLUMNS, KICK_COLUMNS

class BackPopPosteriors():
    def __init__(self, file=None, points=None, log_w=None, log_l=None, var_names=None,
                 blobs=None, var_labels=None):
        """Utility class to handle and analyse posterior samples from BackPop.

        Parameters
        ----------
        file : ``str``, optional
            Path to an HDF5 file containing posterior samples. The file should contain datasets
            named 'points', 'log_w', 'log_l', 'var_names', and 'blobs'.
        points : ``np.ndarray``, optional
            Array of shape (n_samples, n_vars) containing posterior samples.
        log_w : ``np.ndarray``, optional
            Array of shape (n_samples,) containing log weights for each sample.
        log_l : ``np.ndarray``, optional
            Array of shape (n_samples,) containing log likelihoods for each sample.
        var_names : ``list`` of ``str``, optional
            List of variable names corresponding to the columns in `points`.
        blobs : ``np.ndarray``, optional
            Array of shape (n_samples, ...) containing additional data associated with each sample.
            The exact shape and contents depend on the simulation outputs.
        var_labels : ``list`` of ``str``, optional
            List of labels for the variables, used in plotting. If not provided, `var_names` will be used.

        Raises
        ------
        ValueError
            If neither `file` nor all of `points`, `log_w`, `log_l`, and `var_names` are provided.
        """
        # load data from file if provided
        if file is not None:
            with h5.File(file, 'r') as f:
                self.points = f['points'][:]
                self.log_w = f['log_w'][:]
                self.log_l = f['log_l'][:]
                self.var_names = np.array(f['var_names'][:].astype(str).tolist())
                self.blobs = f['blobs'][:]
        # otherwise use provided data
        elif points is not None and log_w is not None and log_l is not None and var_names is not None:
            self.points = points
            self.log_w = log_w
            self.log_l = log_l
            self.var_names = var_names
            self.blobs = blobs 
        # or shout at the user
        else:
            raise ValueError("Must provide either a file or points, log_w, log_l, and labels directly.")

        self.labels = var_labels if var_labels is not None else self.var_names

        # if blobs are provided, parse them into dataframes
        if self.blobs is not None:
            self.bpp = pd.DataFrame(self.blobs["bpp"].reshape(-1, BPP_SHAPE[-1]), columns=BPP_COLUMNS)
            self.kick_info = pd.DataFrame(self.blobs["kick_info"].reshape(-1, KICK_SHAPE[-1]),
                                          columns=KICK_COLUMNS)

            # set index so we can easily filter based on binaries
            self.bpp.index = np.repeat(np.arange(self.bpp.shape[0] / BPP_SHAPE[0]), BPP_SHAPE[0]).astype(int)
            self.kick_info.index = np.repeat(np.arange(self.kick_info.shape[0] / KICK_SHAPE[0]),
                                             KICK_SHAPE[0]).astype(int)

            # filter out empty data (evol_type would never be 0 in a real binary)
            self.bpp = self.bpp[self.bpp["evol_type"] > 0.0]

    def __len__(self):
        return self.points.shape[0]

    def __repr__(self):
        return (f"<BackPopPosteriors: {self.points.shape[0]} samples, "
                f"{self.points.shape[1]} variables")
    
    @property
    def n_vars(self):
        return self.points.shape[1]

    def cornerplot(self, which_vars=None, show=True, **kwargs):
        mask = np.ones(self.n_vars, dtype=bool)
        if which_vars is not None:
            mask = np.isin(self.var_names, which_vars)
            if not np.any(mask):
                raise ValueError("No matching variable names found.")
        fig = corner.corner(
            self.points[:, mask], weights=np.exp(self.log_w), bins=kwargs.pop("bins", 20),
            labels=self.labels[mask], color=kwargs.pop("color", '#074662'),
            plot_datapoints=kwargs.pop("plot_datapoints", False),
            range=kwargs.pop("range", np.repeat(0.999, mask.sum())), **kwargs)

        if show:
            plt.show()
        return fig