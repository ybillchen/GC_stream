from functools import partial

import numpy as np
from scipy.interpolate import RegularGridInterpolator

def gaussian(u):
    return np.exp(-0.5*u**2)

def kde(Xdata, Xeval, h, proj_axes=[], error=None, kernel=gaussian):
    assert kernel == gaussian # currently only support gaussian
    non_proj_axes = np.ones(Xdata.shape[1], dtype=bool)
    non_proj_axes[proj_axes] = False
    if error is None:
        heff = h[:,np.newaxis,non_proj_axes]
    else:
        heff = np.sqrt(
            h[:,np.newaxis,non_proj_axes]**2 + \
            error[np.newaxis,:,non_proj_axes]**2
        )
    dist = np.sqrt(np.sum(
        ((
            Xeval[np.newaxis,:,non_proj_axes] - \
            Xdata[:,np.newaxis,non_proj_axes]
        ) / heff)**2, 
        axis=2
    ))
    k = np.sum(
        kernel(dist) / np.prod(2.50663 * heff,axis=2),  # 2.50663 = sqrt(2*pi)
        axis=0
    )
    return k / len(Xdata)

class PDF(object):
    """docstring for PDF"""
    def __init__(self, data, grids, hs, groups):
        """
        data: array-like (N, M): N data points with M dimensions
        grids: list (M): arrays of evaluation grids
        hs: array-like (N, M): bandwidths for Gaussian KDE
        groups: list
        interp: bool
        """
        self.N, self.dimension = data.shape
        assert self.dimension == len(grids)
        self.data = data
        self.grids = grids
        self.hs = hs
        self.groups = groups

        self.get_pdf()

    def get_pdf(self):
        self.pdfs = []
        for group in self.groups:
            group_data = self.data[:,group]
            group_h = self.hs[:,group]
            group_grid = [self.grids[i] for i in group]
            if group_grid[0] is None:
                # direct KDE
                self.pdfs.append(partial(kde, Xdata=group_data, h=group_h))
            else:
                # interpolation
                group_mesh = np.meshgrid(*group_grid, indexing='ij')
                group_mesh_flatten = np.column_stack([
                    m.flatten() for m in group_mesh
                ])
                
                group_prob = kde(group_data, group_mesh_flatten, group_h)
                group_prob = group_prob.reshape(group_mesh[0].shape)
                self.pdfs.append(RegularGridInterpolator(
                    group_grid,
                    group_prob,
                    bounds_error=False,
                    fill_value=0
                ))

    def eval_pdf(self, data_eval, err_eval=None):
        prob = np.ones(len(data_eval))
        for group, pdf in zip(self.groups, self.pdfs):
            group_grid = [self.grids[i] for i in group]
            if group_grid[0] is None and not err_eval is None:
                prob *= pdf(Xeval=data_eval[:,group], error=err_eval[:,group])
            else:
                prob *= pdf(data_eval[:,group])
        return prob