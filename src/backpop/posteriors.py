import corner
import matplotlib.pyplot as plt
import numpy as np
import h5py as h5

class BackPopPosteriors():
    def __init__(self, file=None, points=None, log_w=None, log_l=None, labels=None, blobs=None):
        
        if file is not None:
            with h5.File(file, 'r') as f:
                self.points = f['points'][:]
                self.log_w = f['log_w'][:]
                self.log_l = f['log_l'][:]
                self.labels = f['labels'][:].astype(str).tolist()
                self.blobs = f['blobs'][:]
        elif points is not None and log_w is not None and log_l is not None and labels is not None:
            self.points = points
            self.log_w = log_w
            self.log_l = log_l
            self.labels = labels
        else:
            raise ValueError("Must provide either a file or points, log_w, log_l, and labels directly.")
        self.blobs = blobs

    def cornerplot(self, show=True, **kwargs):
        fig = corner.corner(
            self.points, weights=np.exp(self.log_w), bins=kwargs.pop("bins", 20),
            labels=self.labels, color=kwargs.pop("color", '#074662'),
            plot_datapoints=kwargs.pop("plot_datapoints", False),
            range=kwargs.pop("range", np.repeat(0.999, len(self.labels))), **kwargs)

        if show:
            plt.show()
        return fig