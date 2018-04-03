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
    """
    See Peter & Rangarajan, "Maximum likelihood wavelet density estimation for image and shape matching" (2008)
    """
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
        self.xs = xs
        self.allocate_vars()
        self.calc_pdf()
        self.calc_grad_nll()
        self.maximise_likelihood()
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

    def allocate_vars(self):
        self.coeffs = {}
        self.nums = {}
        qxs = list(all_qx(self.dim))
        self.var_ix = 0
        self.allocate_vars_j(self.j0, qxs[0:1], 'alpha')
        j0_level_ix = self.var_ix
        for j in range(self.j0, self.j1 + 1):
            self.allocate_vars_j(j, qxs[1:], 'beta')
        self.vars = np.zeros(self.var_ix + 1, dtype=np.float64)

    def allocate_vars_j(self, j, qxs, what):
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
                # store position of coefficient for var_ix
                self.coeffs[j][qx][zs] = self.var_ix
                self.var_ix += 1
                ## print('v%d = \\%s[%d, qx=%s, zs=%s]' % (self.var_ix-1, what, j, str(qx), str(zs)))

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
            return vals * vals
        pdffun.__dict__['dim'] = self.dim
        self.pdf = pdffun

    def calc_grad_nll(self):
        ## self.vars : array with all paramaters, $\alpha_{j,k}$ and $\beta_{j,q,k}$
        ## the Jacobian is the derivative wrt each of those entries.
        ## Jacobian is 2/(g_) * \Sigma \psi_{j,k,q}(x_i), g_ = \sum aprox to \sqrt{f}
        def fun(new_vars, coords):
            jac = np.zeros(len(self.vars))
            for j in self.coeffs:
                jpow2 = 2 ** j
                for qx in self.coeffs[j]:
                    wavef = self.dual_wave_funs[qx]
                    for zs in self.coeffs[j][qx]:
                        ix = self.coeffs[j][qx][zs]
                        jac[ix] = wavef(jpow2, zs, coords).sum()
            g = self.pre_pdf(new_vars, coords).sum()
            return - jac / g
        fun.__dict__['dim'] = self.dim
        self.grad_nll = fun

    def maximise_likelihood(self):

        def fun_grad_nll(new_vars):
            # note this gradient does not include the \frac{\partial h}{\partial x} as per paper
            return self.grad_nll(new_vars, self.xs)

        def fun_nll(new_vars, _lambda):
            vals = self.pre_pdf(new_vars, self.xs)
            return -np.log(vals * vals).sum() + _lambda * fun_h(new_vars)

        def fun_h(new_vars):
            return (new_vars * new_vars).sum() - 1

        def fun_grad_h(new_vars):
            return 2 * new_vars


        print('minimizing', len(self.vars), 'variables and lambda')
        # method='Newton-CG' if we want to do use Lagragian multiplier by hand

        self._lambda = 1
        # random initial vector, subject to constraint
        self.vars = np.random.randn(len(self.vars))
        print(self.vars.shape)
        self.vars = self.vars / ((self.vars * self.vars).sum() ** 0.5)
        max_iters = 1
        nvars = len(self.vars)
        print('N VARS', nvars)
        prev_f = None
        must_run = True
        while must_run and max_iters < 50:

            # see Peter & Rangarajan (2008), eq (27)

            a_vector = fun_grad_h(self.vars)
            a_vector = a_vector.reshape((nvars, 1)) # column vector
            grad_nll = fun_grad_nll(self.vars)
            grad_nll = grad_nll.reshape((nvars, 1)) # column vector
            l_vector = grad_nll + self._lambda * a_vector

            h = fun_h(self.vars)

            # see 2002, Lu, Shiou - Inverses of 2 × 2 Block Matrices, eq (2.2)
            # or http://www.cs.nthu.edu.tw/~jang/book/addenda/matinv/matinv/ (Special case 1)
            k = -0.25 * np.matmul(a_vector.T, a_vector) # Schur complement of H, dim 1x1
            inv_11 = 0.25 * np.eye(a_vector.shape[0]) + 0.0625 / k * np.matmul(a_vector, a_vector.T)
            inv_12 = -0.25 / k * a_vector
            inv_21 = -0.25 / k * a_vector.T
            inv_22 = 1.0 / k

            self_vars = self.vars.reshape((nvars, 1)) - np.matmul(inv_11, l_vector) - inv_12 * h
            self_lambda = self._lambda - np.matmul(inv_21, l_vector) - inv_22 * h

            next_f = fun_nll(self_vars, self_lambda)

            print('(%d)' % max_iters, self._lambda, '->', next_f, h)

            if max_iters <= 2:
                self.vars = self_vars.reshape(-1)
                self._lambda = self_lambda
            else:
                must_run = h > 0.0001
                if must_run:
                    self.vars = self_vars.reshape(-1)
                    self._lambda = self_lambda
            prev_f = next_f
            max_iters += 1

        if h > 0.0001:
            print('Minimisation failed', h)
        else:
            print('Minimised!')
