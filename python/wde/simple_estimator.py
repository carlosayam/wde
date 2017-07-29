from __future__ import division
import pywt
import numpy as np
import itertools as itt
from scipy.interpolate import interp1d
from functools import partial
from .common import *

class SimpleWaveletDensityEstimator(object):
    def __init__(self, wave_name, j0=1, j1=None, thresholding=None):
        self.wave = pywt.Wavelet(wave_name)
        self.j0 = j0
        self.j1 = j1 if j1 is not None else (j0 - 1)
        self.phi_support, self.psi_support = wave_support_info(self.wave)
        self.pdf = None
        if thresholding is None:
            self.thresholding = lambda n, j, dn, c: c
        else:
            self.thresholding = thresholding

    def fit(self, xs):
        "Fit estimator to data. xs is a numpy array of dimension n x d, n = samples, d = dimensions"
        self.dim = xs.shape[1]
        self.dimpow = 2 ** self.dim
        self.calc_wavefuns()
        self.minx = np.amin(xs, axis=0)
        self.maxx = np.amax(xs, axis=0)
        self.n = xs.shape[0]
        self.calc_coefficients(xs)
        self.pdf = self.calc_pdf()
        return True

    def calc_wavefuns(self):
        self.wave_funs = {}
        phi, psi, _ = self.wave.wavefun(level=8)
        phi = interp1d(np.linspace(*self.phi_support, num=len(phi)), phi, fill_value=0.0, bounds_error=False)
        psi = interp1d(np.linspace(*self.psi_support, num=len(psi)), psi, fill_value=0.0, bounds_error=False)
        for wave_x, qx in all_qx(self.dim):
            f = partial(wave_tensor, qx, phi, psi)
            f.qx = qx
            f.support = support_tensor(qx, self.phi_support, self.psi_support)
            f.suppf = partial(suppf_tensor, qx, self.phi_support, self.psi_support)
            self.wave_funs[tuple(qx)] = f

    def calc_coefficients(self, xs):
        self.coeffs = {}
        self.nums = {}
        qxs = list(all_qx(self.dim))
        norm_const = self.do_calculate_j(self.j0, qxs[0:1], xs)
        for j in range(self.j0, self.j1 + 1):
            norm_const += self.do_calculate_j(j, qxs[1:], xs)
        self.norm_const = norm_const

    def do_calculate_j(self, j, qxs, xs):
        jpow2 = 2 ** j
        if j not in self.coeffs:
            self.coeffs[j] = {}
            self.nums[j] = {}
        norm_j = 0.0
        for ix, qx in qxs:
            wavef = self.wave_funs[qx]
            zs_min, zs_max = zs_range(wavef, self.minx, self.maxx, j)
            self.coeffs[j][qx] = {}
            self.nums[j][qx] = {}
            for zs in itt.product(*all_zs_tensor(zs_min, zs_max)):
                v = self.coeffs[j][qx][zs] = calc_coeff_simple(wavef, jpow2, zs, xs)
                self.nums[j][qx][zs] = calc_num(wavef.suppf, jpow2, zs, xs)
                norm_j += v * v
        return norm_j

    def get_betas(self, j):
        return [coeff for ix, qx in list(all_qx(self.dim))[1:] for coeff in self.coeffs[j][qx].values()]

    def get_nums(self):
        return [coeff
                for j in self.nums
                    for ix, qx in list(all_qx(self.dim))[1:]
                        for coeff in self.nums[j][qx].values()]

    def calc_pdf(self):
        def pdffun_j(coords, xs_sum, j, qxs, threshold):
            jpow2 = 2 ** j
            for ix, qx in qxs:
                wavef = self.wave_funs[qx]
                for zs, coeff in self.coeffs[j][qx].iteritems():
                    num = self.nums[j][qx][zs]
                    coeff_t = self.thresholding(self.n, j - self.j0, num, coeff) if threshold else coeff
                    vals = coeff_t * wavef(jpow2, zs, coords)
                    xs_sum += vals
        def pdffun(coords):
            xs_sum = np.zeros(coords[0].shape, dtype=np.float64)
            qxs = list(all_qx(self.dim))
            pdffun_j(coords, xs_sum, self.j0, qxs[0:1], False)
            for j in range(self.j0, self.j1 + 1):
                pdffun_j(coords, xs_sum, j, qxs[1:], True)
            return np.maximum(xs_sum/self.norm_const, 0.0)
        # TODO: this is 2D only
        X = np.linspace(0.0,1.0, num=256)
        Y = np.linspace(0.0,1.0, num=256)
        pred_Z = pdffun(tuple(np.meshgrid(X, Y)))
        full_sum = pred_Z.sum()
        factor = (len(X) * len(Y) / full_sum) if full_sum > 0 else 0.0
        def pdf2(coords):
            return pdffun(coords) * factor
        return pdf2
