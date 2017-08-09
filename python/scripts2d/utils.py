from __future__ import division

import sys
import os
import fcntl
import numpy as np
import scipy.stats as stats
import pandas

from scipy.interpolate import interp2d

class Beta2D(object):
    def __init__(self, a, b, code='beta'):
        self.dist = stats.beta(a, b)
        self.code = code

    def rvs(self, num):
        return self.dist.rvs((num, 2))

    def pdf(self, grid):
        return self.dist.pdf(grid[0]) * self.dist.pdf(grid[1])

class TransformedDirichlet2D(object):
    def __init__(self, alphas, code='diri'):
        if len(alphas) != 3 or min(alphas) <= 0:
            raise ValueError('alphas must be 3 positive numbers')
        self.dist = stats.dirichlet(alpha=alphas)
        self.code = code
        # this code is 2D specific
        z = self._pdf(mise_mesh())
        nns = reduce(lambda x, y: (x-1) * (y-1), z.shape)
        self.sum = z.sum()/nns

    def rvs(self, num):
        x = self.dist.rvs(num)
        return np.stack([np.power(x[:,0],1/4), x[:,1]], axis=1)

    def _pdf(self, grid):
        XX, YY = grid
        XX = np.power(XX, 4)
        ZZ = 1 - XX - YY
        z_pos = (ZZ >= 0) & (ZZ <= 1)
        return self.dist.pdf(ZZ)

    def pdf(self, grid):
        return self._pdf(grid)/self.sum


class TruncatedMultiNormal2D(object):
    def __init__(self, probs, mus, covs, code='mult'):
        self.code = code
        self.probs = probs
        self.dists = [stats.multivariate_normal(mean=mu, cov=cov) for mu, cov in zip(mus, covs)]
        # this code is 2D specific
        z = self._pdf(mise_mesh())
        nns = reduce(lambda x, y: (x-1) * (y-1), z.shape)
        self.sum = z.sum()/nns

    def _rvs(self):
        while True:
            for xvs in zip(*[dist.rvs(100) for dist in self.dists]):
                i = np.random.choice(np.arange(0,len(self.probs)), p=self.probs)
                yield xvs[i]

    def rvs(self, num):
        data = []
        while num > 0:
            for d in self._rvs():
                if 0 <= d[0] and d[0] <= 1 and 0 <= d[1] and d[1] <= 1:
                    data.append(d)
                    num -= 1
                    if num == 0:
                        break
        return np.array(data)

    def _pdf(self, grid):
        pos = np.empty(grid[0].shape + (2,))
        pos[:, :, 0], pos[:, :, 1] = grid
        vals = [dist.pdf(pos) for dist in self.dists]
        pdf_vals = vals[0] * self.probs[0]
        for i in range(len(self.probs) - 1):
            pdf_vals = np.add(pdf_vals, vals[i+1] * self.probs[i+1])
        #pdf_vals = pdf_vals / total
        return pdf_vals

    def pdf(self, grid):
        return self._pdf(grid)/self.sum

def dist_from_code(code):
    if code == 'beta':
        return Beta2D(2, 4, code=code)
    elif code == 'mult':
        sigma = 0.05
        return TruncatedMultiNormal2D(
            [1/9, 8/9],
            [np.array([0.2, 0.3]), np.array([0.7, 0.7])],
            [np.array([[sigma/6, 0], [0, sigma/6]]), np.array([[0.1, sigma/8], [sigma/8, 0.1]])],
            code=code
            )
    elif code == 'mix1':
        sigma = 0.05
        m1 = np.array([[sigma/6, 0], [0, sigma/6.5]])
        return TruncatedMultiNormal2D(
            [1/2, 1/2],
            [np.array([0.2, 0.3]), np.array([0.7, 0.7])],
            [m1, m1],
            code=code
            )
    elif code == 'mix2':
        sigma = 0.05
        angle = 10.
        theta = (angle/180.) * np.pi
        rot = np.array([[np.cos(theta), -np.sin(theta)],
                       [np.sin(theta),  np.cos(theta)]])
        m1 = np.array([[sigma/6, 0], [0, sigma/8]])
        m2 = np.dot(rot, np.dot(m1, rot.T))
        return TruncatedMultiNormal2D(
            [1/2, 1/2],
            [np.array([0.4, 0.3]), np.array([0.7, 0.7])],
            [m1, m2],
            code=code
            )
    elif code == 'mix3':
        sigma = 0.03
        angle = 10.
        theta = (angle/180.) * np.pi
        rot = np.array([[np.cos(theta), -np.sin(theta)],
                       [np.sin(theta),  np.cos(theta)]])
        m1 = np.array([[sigma/6, 0], [0, sigma/7]])
        m2 = np.dot(rot, np.dot(m1, rot.T))
        prop = np.array([8,4,2,1])
        prop = prop/prop.sum()
        return TruncatedMultiNormal2D(
            prop.tolist(),
            [np.array([0.2, 0.3]), np.array([0.5, 0.5]), np.array([0.65, 0.7]), np.array([0.82, 0.85])],
            [m1, m2/2, m1/4, m2/8],
            code=code
            )
    elif code == 'diri':
        return Dirichlet2D([2,3,7])
    else:
        raise NotImplemented('Unknown distribution code [%s]' % code)

def mkdir(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def write_sample(n, i, data):
    mkdir('data2d/samples')
    fname = 'data2d/samples/data-%05d-%03d.csv' % (n, i)
    np.savetxt(fname, data, fmt='%f', delimiter=',')
    return fname

def read_sample(fname):
    return np.genfromtxt(fname, delimiter=',')

def sample_name(fname):
    return os.path.splitext(os.path.basename(fname))[0]

PLAN_FNAME = 'data2d/plan.csv'
def write_plans(plans):
    mkdir('data2d')
    df = pandas.DataFrame(plans)
    print 'sorting & saving, shape=', df.shape
    df = df.sort_values(['rand'])
    df.to_csv(PLAN_FNAME, index=False)

def read_plans(bag_size, bag_number):
    df = pandas.read_csv(PLAN_FNAME)
    start = bag_size * (bag_number - 1)
    rows = df.iloc[start:start + bag_size]
    return rows

def mise_mesh():
    X = np.linspace(0.0,1.0, num=256)
    Y = np.linspace(0.0,1.0, num=256)
    return np.meshgrid(X, Y) # X,Y

def write_dist_pdf(dist):
    mkdir('data2d')
    XX, YY = mise_mesh()
    Z = dist.pdf((XX, YY))
    fname = 'data2d/true-pdf.csv'
    np.savetxt(fname, Z, fmt='%f',delimiter=',')

DIST_PDFS = {}
def read_dist_pdf():
    fname = 'data2d/true-pdf.csv'
    if fname in DIST_PDFS:
        return DIST_PDFS[fname]
    Z = np.genfromtxt(fname, delimiter=',')
    DIST_PDFS[fname] = Z
    return DIST_PDFS[fname]

def calc_ise(pred_pdf, pdf_vals):
    XX, YY = mise_mesh()
    pred_Z = pred_pdf((XX, YY))
    diff = pred_Z - pdf_vals
    # Ok, because we know domain [0,1]x[0,1] => area = 1 x 1 = 1
    err = (diff * diff).sum()
    # extreme values are zero, which do not contribute to integral, hence correction
    # in size "_ - 1".
    nns = reduce(lambda x, y: (x-1) * (y-1), pred_Z.shape)
    return err / nns

def write_ise(fh_ise, fname, dist_code, wave_code, n, j0, j1, k, ise, elapsed_time):
    new_entry = '"%s", "%s", "%s", %d, %d, %d, %d, %f, %f\n' % (fname, dist_code, wave_code, n, j0, j1, k, ise, elapsed_time)
    fh_ise.write(new_entry)

def write_coeffs(fh_coeffs, fname, dist_code, wave_code, n, j0, j1, k, wde):
    for j, coeffs_j in wde.coeffs.iteritems():
        for qx, coeffs_qx in coeffs_j.iteritems():
            for zs, coeff in coeffs_qx.iteritems():
                if abs(coeff) < 0.000001:
                    continue
                new_entry = '"%s", "%s", "%s", %d, %d, %d, %d, %d, "%s", "%s", %f\n' % (fname, dist_code, wave_code, n, j0, j1, k, j, str(qx), str(zs), coeff)
                fh_coeffs.write(new_entry)
