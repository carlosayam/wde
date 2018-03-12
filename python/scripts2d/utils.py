import os
from functools import reduce
import numpy as np
import scipy.stats as stats
import pandas


class Beta2D(object):
    def __init__(self, a, b, code='beta'):
        self.dist = stats.beta(a, b)
        self.code = code
        self.dim = 2

    def rvs(self, num):
        return self.dist.rvs((num, 2))

    def pdf(self, grid):
        return self.dist.pdf(grid[0]) * self.dist.pdf(grid[1])

class TruncatedMultiNormalD(object):
    """Truncated mixture or multivariate normal distributions. Dimension is inferred from first $\mu$"""

    def __init__(self, probs, mus, covs, code='mult'):
        self.code = code
        self.probs = probs
        self.dists = [stats.multivariate_normal(mean=mu, cov=cov) for mu, cov in zip(mus, covs)]
        self.dim = len(mus[0])
        z = self._pdf(mise_mesh(self.dim))
        nns = reduce(lambda x, y: (x-1) * (y-1), z.shape)
        self.sum = z.sum()/nns

    def mathematica(self):
        # render Mathematica code to plot
        def fn(norm_dist):
            mu = np.array2string(norm_dist.mean, separator=',')
            mu = mu.replace('[','{').replace(']','}').replace('e','*^')
            cov = np.array2string(norm_dist.cov, separator=',')
            cov = cov.replace('[','{').replace(']','}').replace('e','*^')
            return 'MultinormalDistribution[%s,%s]' % (mu, cov)
        probs = '{%s}' % ','.join([str(f / min(self.probs)) for f in self.probs])
        dists = '{%s}' % ','.join([fn(d) for d in self.dists])
        resp = 'MixtureDistribution[%s,%s]' % (probs, dists)
        return resp

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
        if type(grid) == tuple or type(grid) == list:
            pos = np.stack(grid, axis=0)
            pos = np.moveaxis(pos, 0, -1)
        else:
            pos = grid
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
    elif code == 'mult' or code == 'mul2':
        sigma = 0.01
        return TruncatedMultiNormalD(
            [1.5/9, 7.5/9],
            [np.array([0.2, 0.3]), np.array([0.7, 0.7])],
            [np.array([[sigma/6, 0], [0, sigma/6]]), np.array([[0.015, sigma/64], [sigma/64, 0.015]])],
            code=code
        )
    elif code == 'mul3':
        sigma = 0.01
        return TruncatedMultiNormalD(
            [0.4, 0.3, 0.3],
            [np.array([0.3, 0.4, 0.35]),
             np.array([0.7, 0.7, 0.6]),
             np.array([0.7, 0.6, 0.35])],
            [np.array([[0.02, 0.01, 0.], [0.01, 0.02, 0.], [0., 0., 0.02]]),
             np.array([[0.0133333, 0., 0.], [0., 0.0133333, 0.], [0., 0., 0.0133333]]),
             np.array([[0.025, 0., 0.], [0., 0.025, 0.01], [0., 0.01, 0.025]])
             ],
            code=code
            )
    elif code == 'mix1':
        sigma = 0.05
        m1 = np.array([[sigma/6, 0], [0, sigma/6.5]])
        return TruncatedMultiNormalD(
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
        return TruncatedMultiNormalD(
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
        return TruncatedMultiNormalD(
            prop.tolist(),
            [np.array([0.2, 0.3]), np.array([0.5, 0.5]), np.array([0.65, 0.7]), np.array([0.82, 0.85])],
            [m1, m2/2, m1/4, m2/8],
            code=code
            )
    elif code == 'mix4':
        sigma = 0.03
        angle = 10.
        theta = (angle / 180.) * np.pi
        rot = np.array([[np.cos(theta), -np.sin(theta)],
                        [np.sin(theta), np.cos(theta)]])
        m1 = np.array([[sigma / 6, 0], [0, sigma / 7]])
        m2 = np.dot(rot, np.dot(m1, rot.T))
        prop = np.array([8, 4, 2, 1, 384])
        prop = prop / prop.sum()
        return TruncatedMultiNormalD(
            prop.tolist(),
            [np.array([0.2, 0.3]), np.array([0.5, 0.5]), np.array([0.65, 0.7]), np.array([0.82, 0.85]), np.array([0.5, 0.5])],
            [m1, m2 / 2, m1 / 4, m2 / 8, 0.18 * np.eye(2, 2)],
            code=code
        )
    elif code == 'mix5':
        sigma = 0.03
        angle = 10.
        theta = (angle / 180.) * np.pi
        rot = np.array([[np.cos(theta), -np.sin(theta)],
                        [np.sin(theta), np.cos(theta)]])
        m1 = np.array([[sigma / 6, 0], [0, sigma / 7]])
        m2 = np.dot(rot, np.dot(m1, rot.T))
        prop = np.array([8, 4, 2, 1])
        prop = prop / prop.sum()
        return TruncatedMultiNormalD(
            prop.tolist(),
            [np.array([0.2, 0.3]), np.array([0.5, 0.5]), np.array([0.65, 0.7]), np.array([0.82, 0.85]),
             np.array([0.5, 0.5])],
            [m1, m2 / 2, m1 / 6, m2 / 8],
            code=code
        )
    elif code == 'mix6':
        theta = np.pi / 4
        rot = lambda angle : np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        m0 = np.array([[0.1, 0], [0, 0.0025]])
        m1 = np.dot(rot(theta), np.dot(m0, rot(theta).T)) / 2
        m2 = np.dot(rot(-theta), np.dot(m0, rot(-theta).T)) / 2
        prop = np.array([1, 1])
        prop = prop / prop.sum()
        return TruncatedMultiNormalD(
            prop.tolist(),
            [np.array([0.3, 0.3]), np.array([0.7, 0.3])], [m1, m2],
            code=code
        )
    elif code == 'mix7': ## not good
        m0 = np.array([[0.1, 0], [0, 0.005]])
        m1 = np.array([[0.005, 0], [0, 0.1]])
        prop = np.array([1, 1, 1, 1])
        prop = prop / prop.sum()
        return TruncatedMultiNormalD(
            prop.tolist(),
            [np.array([0.5, 0.3]),
             np.array([0.5, 0.7]),
             np.array([0.3, 0.5]),
             np.array([0.7, 0.5])],
            [m0, m0, m1, m1],
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
    ## np.savetxt(fname, data, fmt='%f', delimiter=',') !!!
    return fname

def read_sample(fname):
    return np.genfromtxt(fname, delimiter=',')

def sample_name(fname):
    return os.path.splitext(os.path.basename(fname))[0]

PLAN_FNAME = 'data2d/plan.csv'
def write_plans(plans):
    mkdir('data2d')
    df = pandas.DataFrame(plans)
    print('sorting & saving, shape=', df.shape)
    df = df.sort_values(['rand'])
    df.to_csv(PLAN_FNAME, index=False)

def read_plans(bag_size, bag_number):
    df = pandas.read_csv(PLAN_FNAME)
    start = bag_size * (bag_number - 1)
    rows = df.iloc[start:start + bag_size]
    return rows

def mise_mesh(d=2):
    grid_n = 256 if d == 2 else 64
    VVs = [np.linspace(0.0,1.0, num=grid_n) for i in range(d)]
    return np.meshgrid(*VVs)

def write_dist_pdf(dist):
    mkdir('data2d')
    Z = dist.pdf(mise_mesh(dist.dim))
    fname = 'data2d/true-pdf.npy'
    np.save(fname, Z, allow_pickle=False)

DIST_PDFS = {}
def read_dist_pdf():
    global DIST_PDFS
    fname = 'data2d/true-pdf.npy'
    if fname in DIST_PDFS:
        return DIST_PDFS[fname]
    Z = np.load(fname, allow_pickle=False)
    DIST_PDFS[fname] = Z
    return DIST_PDFS[fname]

def calc_ise(pred_pdf, pdf_vals):
    pred_Z = pred_pdf(mise_mesh(d=pred_pdf.dim))
    diff = pred_Z - pdf_vals
    # Ok, because we know domain [0,1]x[0,1] => area = 1 x 1 = 1
    err = (diff * diff).sum()
    # extreme values are zero, which do not contribute to integral, hence correction
    # in size "_ - 1".
    nns = reduce(lambda x, y: (x - 1) * (y - 1), pred_Z.shape)
    return err / nns

def calc_hellinger(pred_pdf, pdf_vals):
    pred_Z = np.sqrt(np.clip(pred_pdf(mise_mesh(d=pred_pdf.dim)), a_min=0, a_max=None))
    diff = pred_Z - np.sqrt(pdf_vals)
    # Ok, because we know domain [0,1]x[0,1] => area = 1 x 1 = 1
    err = (diff * diff).sum()
    # extreme values are zero, which do not contribute to integral, hence correction
    # in size "_ - 1".
    nns = reduce(lambda x, y: (x-1) * (y-1), pred_Z.shape)
    return err / nns

def write_ise(fh_ise, fname, dist_code, wave_code, n, j0, j1, k, ise, hd, elapsed_time):
    new_entry = '"%s", "%s", "%s", %d, %d, %d, %d, %f, %f, %f\n' % (fname, dist_code, wave_code, n, j0, j1, k, ise, hd, elapsed_time)
    print(new_entry)
    fh_ise.write(new_entry)

def write_coeffs(fh_coeffs, fname, dist_code, wave_code, n, j0, j1, k, wde):
    for j, coeffs_j in wde.coeffs.items():
        for qx, coeffs_qx in coeffs_j.items():
            for zs, coeff in coeffs_qx.items():
                if abs(coeff) < 0.000001:
                    continue
                new_entry = '"%s", "%s", "%s", %d, %d, %d, %d, %d, "%s", "%s", %f\n' % (fname, dist_code, wave_code, n, j0, j1, k, j, str(qx), str(zs), coeff)
                fh_coeffs.write(new_entry)
