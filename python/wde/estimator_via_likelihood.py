import pywt
import numpy as np
import math
import itertools as itt
from scipy.interpolate import interp1d
from functools import partial
from .common import wave_support_info, all_qx, wave_tensor, support_tensor, suppf_tensor, calculate_nearest_balls, \
    zs_range, all_zs_tensor, calc_num, calc_coeff
from scipy.optimize import minimize

class WaveletDensityEstimatorByLikelihood(object):
    def __init__(self, wave_name, k=1, j0=1, j1=None):
        self.wave = pywt.Wavelet(wave_name)
        self.k = k
        self.j0 = j0
        self.j1 = j1 if j1 is not None else (j0 - 1)
        self.multi_supports = wave_support_info(self.wave)
        self.pdf = None

    def fit(self, xs):
        "Fit estimator to data. xs is a numpy array of dimension n x d, n = samples, d = dimensions"
        self.dim = xs.shape[1]
        self.dimpow = 2 ** self.dim
        self.set_wavefuns(self.dim)
        self.minx = np.amin(xs, axis=0)
        self.maxx = np.amax(xs, axis=0)
        self.n = xs.shape[0]
        self.allocate_vars(xs)
        self.calc_pdf()
        self.calc_jac()
        self.maximise_likelihood(xs)
        return True

    def set_wavefuns(self, dim):
        self.wave_funs = self.calc_wavefuns(dim, 'base', self.multi_supports, self.wave)
        self.dual_wave_funs = self.calc_wavefuns(dim, 'dual', self.multi_supports, self.wave)

    @staticmethod
    def calc_wavefuns(dim, which, supports, wave):
        resp = {}
        phi_support, psi_support = supports[which]
        defs = wave.wavefun(level=12)
        if len(defs) == 5:
            if which == 'base':
                phi, psi = defs[0], defs[1]
            else:
                phi, psi = defs[2], defs[3]
        else:
            phi, psi = defs[0], defs[1]
        phi = interp1d(np.linspace(*phi_support, num=len(phi)), phi, fill_value=0.0, bounds_error=False)
        psi = interp1d(np.linspace(*psi_support, num=len(psi)), psi, fill_value=0.0, bounds_error=False)
        for wave_x, qx in all_qx(dim):
            f = partial(wave_tensor, qx, phi, psi)
            f.qx = qx
            f.support = support_tensor(qx, phi_support, psi_support)
            f.suppf = partial(suppf_tensor, qx, phi_support, psi_support)
            resp[tuple(qx)] = f
        return resp

    def allocate_vars(self, xs):
        self.coeffs = {}
        self.nums = {}
        qxs = list(all_qx(self.dim))
        self.var_ix = 0
        self.allocate_vars_j(self.j0, qxs[0:1], xs)
        j0_level_ix = self.var_ix
        for j in range(self.j0, self.j1 + 1):
            self.allocate_vars_j(j, qxs[1:], xs)
        self.vars = np.zeros(self.var_ix + 1, dtype=np.float64)
        # set first level to 1s, leave all other levels to 0 and normalise square to 1
        self.vars[:j0_level_ix] = 1.0
        self.vars = self.vars / math.sqrt(j0_level_ix)
        self.xs = xs

    def allocate_vars_j(self, j, qxs, xs):
        """
        Allocate coefficients for level J into a position in vars array by simply going through everyone of them
        and setting the current index for that coefficient in the vars array
        :param j: level for which coefficients are allocated
        :param qxs: tensor product
        :param xs: values so we can calculate support range
        """
        if j not in self.coeffs:
            self.coeffs[j] = {}
        for ix, qx in qxs:
            wavef = self.wave_funs[qx]
            zs_min, zs_max = zs_range(wavef, self.minx, self.maxx, j)
            self.coeffs[j][qx] = {}
            for zs in itt.product(*all_zs_tensor(zs_min, zs_max)):
                self.coeffs[j][qx][zs] = self.var_ix
                self.var_ix += 1

    def calc_pdf(self):
        def pdffun_j(vars, coords, xs_sum, j, qxs, threshold):
            jpow2 = 2 ** j
            for ix, qx in qxs:
                wavef = self.dual_wave_funs[qx]
                for zs, var_ix in self.coeffs[j][qx].items():
                    coeff = vars[var_ix]
                    vals = coeff * wavef(jpow2, zs, coords)
                    xs_sum += vals
        def pre_pdffun(vars, coords):
            if type(coords) == tuple or type(coords) == list:
                xs_sum = np.zeros(coords[0].shape, dtype=np.float64)
            else:
                xs_sum = np.zeros(coords.shape[0], dtype=np.float64)
            qxs = list(all_qx(self.dim))
            pdffun_j(vars, coords, xs_sum, self.j0, qxs[0:1], False)
            for j in range(self.j0, self.j1 + 1):
                pdffun_j(vars, coords, xs_sum, j, qxs[1:], True)
            ## Note: norm_const = self.vars * self.vars -- should be close to 1 by constraint
            return xs_sum
        self.pre_pdf = pre_pdffun
        def pdffun(coords):
            vals = pre_pdffun(self.vars, coords)
            return (vals * vals).sum()
        pdffun.__dict__['dim'] = self.dim
        self.pdf = pdffun

    def calc_jac(self, new_vars):
        ## self.vars : array with all paramaters, $\alpha_{j,k}$ and $\beta_{j,q,k}$
        ## the Jacobian is the derivative wrt each of those entries.
        ## Jacobian is 2/(g_) * \Sigma \psi_{j,k,q}(x_i), g_ = aprox to \sqrt{f}
        def jacfun(coords):
            if type(coords) == tuple or type(coords) == list:
                xs_sum = np.zeros(coords[0].shape, dtype=np.float64)
            else:
                xs_sum = np.zeros(coords.shape[0], dtype=np.float64)
            qxs = list(all_qx(self.dim))
            pdffun_j(coords, xs_sum, self.j0, qxs[0:1], False)
            for j in range(self.j0, self.j1 + 1):
                pdffun_j(coords, xs_sum, j, qxs[1:], True)
            ## Note: norm_const = self.vars * self.vars -- should be close to 1 by constraint
            return (xs_sum * xs_sum)
        jacfun.__dict__['dim'] = self.dim
        return jacfun

    def maximise_likelihood(self, data_xs):
        id_matrix = np.identity(len(self.vars)) * 4
        def jacfun(new_vars):
            return self.jacobian(new_vars, self.xs)
        def fun(new_vars):
            vals = self.pre_pdf(new_vars, self.xs)
            return -np.log(vals * vals).mean()
        def fun_c(new_vars):
            return (new_vars * new_vars).sum() - 1
        constraints = [
            {
                'type': 'eq',
                'fun': fun_c
            }
        ]
        hessian = lambda new_vars: id_matrix
        print('minimizing', len(self.vars), 'variables')
        resp = minimize(fun, self.vars.copy(), method='SLSQP', jac=jacfun, hess=hessian, constraints=constraints)
        if resp is None or not resp.success:
            raise ValueError('Minimisation failed')
