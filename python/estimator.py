import pywt
import math
import numpy as np
import itertools as itt
from sklearn.neighbors import BallTree
from scipy.special import gamma

class WaveletDensityEstimator(object):
    def __init__(self, wave_name, k = 1, j0 = 0, j1 = 0):
        self.wave = pywt.Wavelet(wave_name)
        self.k = k
        self.j0 = j0
        self.j1 = j1

    def fit(self, xs):
        "Fit estimator to data. xs is a numpy array of dimension n x d, n = samples, d = dimensions"
        self.dim = xs.shape[1]
        self.dimpow = 2 ** self.dim
        self.build_coeff_grid(xs)
        self.calc_coefficients(xs)
        self.fun_def()
        return True

    def build_coeff_grid(self, xs):
        self.minx = np.amin(xs, axis=0)
        self.maxx = np.amax(xs, axis=0)
        self.ks = np.zeros((self.j1 - self.j0 + 1, self.dim, 2), dtype=np.int)
        for j in range(self.j0, self.j1 + 1):
            for d in range(self.dim):
                self.ks[j - self.j0, d,:] = calc_k_min_max(self.wave, j, self.minx[d], self.maxx[d])
                #print 'ks @',j,'d=' + str(d),'=',self.ks[j - self.j0, d,:]
        #print 'ks shape', self.ks.shape
        self.coeffs = {}
        for j in range(self.j0, self.j1 + 1):
            shape = (self.dimpow,) + tuple(self.ks[j - self.j0,:,1] - self.ks[j - self.j0,:,0])
            self.coeffs[j] = np.zeros(shape)
            #print 'coeffs shape @', j, self.coeffs[j].shape

    def calc_coefficients(self, xs):
        grid_xs = self.gridify_xs(xs)
        xs_and_nearest = self.calculate_nearest(xs)
        self.coeffs = self.do_calculate(grid_xs, xs_and_nearest)

    def gridify_xs(self, xs):
        grid_xs = {}
        for j in range(self.j0, self.j1+1):
            jpow = 2 ** j
            grid_xs[j] = {}
            if j == self.j0:
                iters = [xrange(int(jpow * self.minx[d]), int(jpow * self.maxx[d]) + 2) for d in range(self.dim)]
                for ks in itt.product(*iters):
                    cond = (np.floor(jpow * xs) == ks).all(axis=1)
                    grid_xs[j][ks] = np.where(cond)
            else:
                for ks_up, where_xs in grid_xs[j-1].iteritems():
                    # TODO theory - one could stop splitting for len <= N0, what does this mean?
                    if len(where_xs[0]) == 0:
                        continue
                    sub_xs = xs[where_xs]
                    ks_up_arr = np.array(ks_up)
                    for ix2s in itt.product(*itt.tee(xrange(2), self.dim)):
                        ks = 2 * ks_up_arr + np.array(ix2s)
                        cond = (np.floor(sub_xs * jpow) == ks).all(axis=1)
                        grid_xs[j][tuple(ks.tolist())] = (where_xs[0][cond],)
        return grid_xs

    def calculate_nearest(self, xs):
        ball_tree = BallTree(xs)
        k_near_balls = ball_tree.query(xs, self.k + 1)[0][:,[-1]]
        factor = calc_factor(xs.shape[0], self.dim, self.k)
        xs[:,:-1] = np.power(k_near_balls, self.dim/2.0) * factor
        return xs

    def do_calculate(self, grid_xs, xs_and_nearest):
        self.coeffs = {}
        for j in range(self.j0, self.j1 + 1):
            shape = (self.dimpow,) + tuple(self.ks[j - self.j0,:,1] - self.ks[j - self.j0,:,0])
            self.coeffs[j] = np.zeros(shape)
            wave_range = xrange(1, self.dimpow) if j == self.j0 else xrange(0,1)
            for wave_x in wave_range:
                pass
        pass

    def fun_def(self):
        pass

def calc_k_min_max(wave, j, minx, maxx):
    jp = 2 ** j
    return [jp * minx - wave.dec_len, jp * maxx + wave.dec_len]

def calc_factor(l, dim, k):
    v_unit = np.pi ** (dim/2) / gamma(dim/2.0 + 1)
    return math.sqrt(v_unit) * gamma(k) / gamma(k + 0.5) / math.sqrt(l)