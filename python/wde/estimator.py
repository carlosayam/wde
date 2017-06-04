from __future__ import division
import pywt
import math
import numpy as np
import itertools as itt
from sklearn.neighbors import BallTree
from scipy.interpolate import interp1d
from scipy.special import gamma
from functools import partial

def np_mult(cols):
    if len(cols) == 1:
        return cols[0]
    if len(cols) == 2:
        return np.multiply(*cols)
    else:
        return np.multiply(cols[0], np_mult(cols[1:]))

# all possible combinations of {0,1}, i.e. {0,1}^dim
def all_qx(dim):
    for wave_x, qx in enumerate(itt.product(xrange(2), repeat=dim)):
        yield wave_x, tuple(qx)

def support_tensor(qx, phi_supp, psi_supp):
    return np.array([
        [(phi_supp if q == 0 else psi_supp)[0] for d,q in enumerate(qx)],
        [(phi_supp if q == 0 else psi_supp)[1] for d,q in enumerate(qx)]
    ])
# all z-min, z-max ranges as tensor of ranges for phi and psi support
def all_zs_tensor(zs_min, zs_max):
    it = zip(zs_min, zs_max)
    return [xrange(int(a), int(b)+1) for a, b in it]

# tensor product of z-min per dimension
def z0_tensor(qx, zs_phi, zs_psi):
    return [(zs_phi if q2 == 0 else zs_psi)[0][d] for d, q2 in enumerate(qx)]

# tensor product of max k per dimension
def z1_tensor(qx, zs_phi, zs_psi):
    return [(zs_phi if q2 == 0 else zs_psi)[1][d] for d, q2 in enumerate(qx)]
    
# tensor product waves
# xs = rows x dim
# zs = dim
def wave_tensor(qx, phi, psi, jpow, zs, xs):
    cols = []
    if type(xs) == tuple:
        proj = lambda xs,i: xs[i]
    else:
        proj = lambda xs,i: xs[:,i]
    for i,q2 in enumerate(qx):
        wavef = phi if q2 == 0 else psi
        xs_proj = proj(xs,i)
        cols.append(wavef(jpow * xs_proj - zs[i]))
    return np_mult(tuple(cols)) * (jpow ** (len(qx)/2.0))

def suppf_tensor(qx, phi_sup, psi_sup, jpow, zs, xs):
    cols = []
    if type(xs) == tuple:
        proj = lambda xs,i: xs[i]
    else:
        proj = lambda xs,i: xs[:,i]
    for i,q2 in enumerate(qx):
        supp = phi_sup if q2 == 0 else psi_sup
        xs_proj = proj(xs,i)
        xjz = jpow * xs_proj - zs[i]
        cols.append(np.greater_equal(xjz, supp[0]))
        cols.append(np.less_equal(xjz, supp[1]))
    return np_mult(tuple(cols))


# factor for num samples l, dimension dim and nearest index k
def calc_factor(l, dim, k):
    v_unit = (np.pi ** (dim/2.0)) / gamma(dim/2.0 + 1)
    return math.sqrt(v_unit) * (gamma(k) / gamma(k + 0.5)) / math.sqrt(l)

# calculate V(k);i for each row xs[i] and return dataset with that attached
def calculate_nearest_balls(k, xs):
    dim = xs.shape[1]
    ball_tree = BallTree(xs)
    k_near_radious = ball_tree.query(xs, k + 1)[0][:,[-1]]
    factor = calc_factor(xs.shape[0], dim, k)
    return np.power(k_near_radious, dim/2.0) * factor

def wave_support_info(wave):
    if wave.family_name in ['Daubechies', 'Symlets']:
        phi_support = (0, wave.dec_len - 1)
        psi_support = (1 - wave.dec_len // 2, wave.dec_len // 2)
    elif wave.family_name in ['Coiflets']:
        phi_support = (1 - wave.dec_len // 2, wave.dec_len // 2)
        psi_support = (1 - wave.dec_len // 2, wave.dec_len // 2)
    else:
        raise ValueError('wave family %s not known support' % wave.family_name)
    return phi_support, psi_support

def gridify_xs(j0, j1, xs, minx, maxx):
    grid_xs = {}
    dim = xs.shape[1]
    for j in range(j0, j1+1):
        jpow = 2 ** j
        grid_xs[j] = {}
        if j == j0:
            iters = [xrange(int(jpow * minx[d]), int(jpow * maxx[d]) + 1) for d in range(dim)]
            for zs in itt.product(*iters):
                cond = (np.floor(jpow * xs) == zs).all(axis=1)
                grid_xs[j][zs] = np.where(cond)
        else:
            for zs_up, where_xs in grid_xs[j-1].iteritems():
                # TODO theory - one could stop splitting for len <= N0, what does this mean?
                if len(where_xs[0]) == 0:
                    continue
                sub_xs = xs[where_xs]
                zs_up_arr = np.array(zs_up)
                for _, ix2s in all_qx(dim):
                    zs = 2 * zs_up_arr + np.array(ix2s)
                    cond = (np.floor(sub_xs * jpow) == zs).all(axis=1)
                    grid_xs[j][tuple(zs.tolist())] = (where_xs[0][cond],)
    return grid_xs

def zs_range(wavef, minx, maxx, j):
    zs_min = np.ceil((2 ** j) * minx - wavef.support[1])
    zs_max = np.floor((2 ** j) * maxx - wavef.support[0])
    return zs_min, zs_max
    
def calc_coeff(wave_tensor_qx, jpow, zs, xs, xs_balls):
    vals = wave_tensor_qx(jpow, zs, xs)
    all_prods = vals * xs_balls[:,0]
    return all_prods.sum()

def soft_threshold(threshold):
    def f(j, n, coeff):
        lvl_factor = math.sqrt(j / n)
        lvl_t = threshold * lvl_factor
        if coeff < 0:
            if -coeff < lvl_t:
                return 0
            else:
                return coeff + lvl_t
        else:
            if coeff < lvl_t:
                return 0
            else:
                return coeff - lvl_t
    return f

def hard_threshold(threshold):
    def f(j, n, coeff):
        lvl_factor = math.sqrt(j / n)
        lvl_t = threshold * lvl_factor
        if coeff < 0:
            if -coeff < lvl_t:
                return 0
            else:
                return coeff
        else:
            if coeff < lvl_t:
                return 0
            else:
                return coeff
    return f

class WaveletDensityEstimator(object):
    def __init__(self, wave_name, k=1, j0=1, j1=None, thresholding=None):
        self.wave = pywt.Wavelet(wave_name)
        self.k = k
        self.j0 = j0
        self.j1 = j1 if j1 is not None else (j0 - 1)
        self.phi_support, self.psi_support = wave_support_info(self.wave)
        if thresholding is None:
            self.thresholding = lambda n, j, c: c
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
        #grid_xs = gridify_xs(self.j0, self.j1, xs, self.minx, self.maxx)
        xs_balls = calculate_nearest_balls(self.k, xs)
        self.do_calculate(xs, xs_balls)

    def do_calculate(self, xs, xs_balls):
        self.coeffs = {}
        qxs = list(all_qx(self.dim))
        norm_const = self.do_calculate_j(self.j0, qxs[0:1], xs, xs_balls)
        for j in range(self.j0, self.j1 + 1):
            norm_const += self.do_calculate_j(j, qxs[1:], xs, xs_balls)
        self.norm_const = norm_const

    def do_calculate_j(self, j, qxs, xs, xs_balls):
        jpow2 = 2 ** j
        if j not in self.coeffs:
            self.coeffs[j] = {}
        norm_j = 0.0
        for ix, qx in qxs:
            wavef = self.wave_funs[qx]
            zs_min, zs_max = zs_range(wavef, self.minx, self.maxx, j)
            self.coeffs[j][qx] = {}
            for zs in itt.product(*all_zs_tensor(zs_min, zs_max)):
                v = self.coeffs[j][qx][zs] = calc_coeff(wavef, jpow2, zs, xs, xs_balls)
                norm_j += v * v
        return norm_j

    def get_betas(self, j):
        return [coeff for ix, qx in list(all_qx(self.dim))[1:] for coeff in self.coeffs[j][qx].values()]

    def calc_pdf(self):
        def pdffun_j(coords, xs_sum, j, qxs, threshold):
            jpow2 = 2 ** j
            for ix, qx in qxs:
                wavef = self.wave_funs[qx]
                for zs, coeff in self.coeffs[j][qx].iteritems():
                    coeff_t = self.thresholding(self.n, j, coeff) if threshold else coeff
                    vals = coeff_t * wavef(jpow2, zs, coords)
                    xs_sum += vals
        def pdffun(coords):
            xs_sum = np.zeros(coords[0].shape, dtype=np.float64)
            qxs = list(all_qx(self.dim))
            pdffun_j(coords, xs_sum, self.j0, qxs[0:1], False)
            for j in range(self.j0, self.j1 + 1):
                pdffun_j(coords, xs_sum, j, qxs[1:], True)
            return (xs_sum * xs_sum)/self.norm_const
        X = np.linspace(0.0,1.0, num=75)
        Y = np.linspace(0.0,1.0, num=75)
        pred_Z = pdffun(tuple(np.meshgrid(X, Y)))
        factor = (len(X) * len(Y) / pred_Z.sum())
        def pdf2(coords):
            return pdffun(coords) * factor
        return pdf2
