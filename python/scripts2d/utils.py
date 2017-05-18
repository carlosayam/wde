class Beta2D(object):
    def __init__(self, code='beta', *args):
        self.args = args
        self.code = code

    def rvs(self, num):
        return np.random.beta(*self.args, (num, 2))

    def pdf(self, grid):
        #print X.shape, X
        #print Y.shape, Y
        return statsm.beta(*self.args, grid)

class TruncatedMultiNormal2D(object):
    def __init__(self, probs, mus, covs, code=code):
        self.code = code
        self._pdf = None
        self.probs = probs
        self.mus = mus
        self.covs = covs
        self.dists = [stats.multivariate_normal(mean=mu, cov=cov) for mu, cov in zip(self.mus, self.covs)]

    def rvs(self, num):
        data = []
        while num > 0:
            for d in mix.mv_mixture_rvs(self.probs, 100, self.dists, 2):
                if 0 <= d[0] and d[0] <= 1 and 0 <= d[1] and d[1] <= 1:
                    data.append(d)
                    num -= 1
                    if num == 0: break
        return np.array(data)

    def pdf(self, grid):
        if self._pdf is None:
            self._pdf = self.calc_pdf()
        return self._pdf(grid)

    def calc_pdf(self):
        mesh = mise_mesh()
        vals = [dist.pdf(mesh) for dist in self.dists]
        pdf_vals = np.add(*[val[i] * self.probs[i] for i in range(len(self.probs))])
        total = pdf_vals.sum()
        pdf_vals = pdf_vals / total
        return interp2d(pdf_vals, fill=0.0)

def dist_from_code(code):
    if code == 'beta':
        return Beta2D(2, 4, code=code)
    elif code == 'mult':
        return TruncatedMultiNormal2D(
            [1/9, 8/9],
            [np.array([0.2, 0.3]), np.array([0.7, 0.7])],
            [np.array([[sigma/6, 0], [0, sigma/6]]), np.array([[0.1, sigma/8], [sigma/8, 0.1]])],
            code=code
            )
    else:
        raise NotImplemented('Unknown distribution code [%s]' % code)

def write_sample(dist, n, i, data):
    fname = 'data2d/%s/samples/data-%05d-%03d.csv' % (dist.code, n, i)
    np.savetxt(fname, data, fmt='%f', delimiter=',')
    return fname

def read_sample(fname):
    data = np.readtxt(fname)
    code = fname.split('/')[1]
    return data, dist_from_code(code)

PLAN_FNAME = 'data2d/plan.csv'
def write_plans(plans):
    df = pandas.DataFrame(plans)
    df = df.sort_values(['j0','j1','k','wave_name','fname'])
    df.to_csv(PLAN_FNAME, index=False)

def read_plans(bag_size, bag_number):
    df = pandas.from_csv(PLAN_FNAME)
    start = bag_size * (bag_number - 1)
    rows = df[start:start + bag_size]
    return rows

def mise_mesh():
    X = np.linspace(0.0,1.0, num=75)
    Y = np.linspace(0.0,1.0, num=75)
    return np.meshgrid(X, Y) # X,Y

def write_dist_pdf(dist):
    X, Y = mise_mesh()
    Z = dist.pdf((X, Y))
    fname = 'data2d/%s/pdf.csv' % dist.code
    np.savetxt(fname, Z, fmt='%f', delimiter=',')

def calc_mise(true):
    pass
