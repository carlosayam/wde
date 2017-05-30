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

class TruncatedMultiNormal2D(object):
    def __init__(self, probs, mus, covs, code='mult'):
        self.code = code
        self._pdf = None
        self.probs = probs
        self.dists = [stats.multivariate_normal(mean=mu, cov=cov) for mu, cov in zip(mus, covs)]

    def _rvs(self):
        while True:
            i = np.random.choice(np.arange(0,len(self.probs)), p=self.probs)
            dist = self.dists[i]
            yield dist.rvs(1)

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

    def pdf(self, grid):
        pos = np.empty(grid[0].shape + (2,))
        pos[:, :, 0], pos[:, :, 1] = grid
        vals = [dist.pdf(pos) for dist in self.dists]
        pdf_vals = np.add(*[vals[i] * self.probs[i] for i in range(len(self.probs))])
        total = pdf_vals.sum()
        #pdf_vals = pdf_vals / total
        return pdf_vals

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
    X = np.linspace(0.0,1.0, num=75)
    Y = np.linspace(0.0,1.0, num=75)
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
    X = np.linspace(0.0,1.0, num=75)
    Y = np.linspace(0.0,1.0, num=75)
    pred_Z = pred_pdf(tuple(np.meshgrid(X, Y)))
    diff = pred_Z - pdf_vals
    err = (diff * diff).sum()
    return err / (len(X) * len(Y))

def write_ise(fhandle, fname, dist_code, wave_code, n, j0, j1, k, ise, elapsed_time):
    new_entry = '"%s", "%s", "%s", %d, %d, %d, %d, %f, %f\n' % (fname, dist_code, wave_code, n, j0, j1, k, ise, elapsed_time)
    fhandle.write(new_entry)
