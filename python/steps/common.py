import os
import numpy as np
import math

if 'PBS_O_HOME' in os.environ:
    PBS_O_HOME = os.environ['PBS_O_HOME']
else:
    PBS_O_HOME = 'data'

class Adder(object):
    def __init__(self, what):
        self.what = what
        self.sum = 0
        self.sum2 = 0
        self.n = 0
    def add_v(self, v):
        self.sum += v
        self.sum2 += v*v
        self.n += 1
    def __str__(self):
        if self.n == 0:
            return '%s: undef' % self.what
        mean = self.sum/self.n
        var = self.sum2/self.n - mean * mean
        return '%s: aver %f (+/- %f)' % \
               (self.what, mean, 2 * math.sqrt(var))


def ensure_dir(directory, overwrite=True):
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except FileExistsError as ex:
            # ignore race condition between multiple jobs
            pass
    else:
        if not overwrite:
            raise ValueError('Directory %s already exists' % directory)

def parent_dir(code):
    return '%(home)s/STEPS/%(dist_code)s' % dict(home=PBS_O_HOME, dist_code=code)

def true_pdf_fname(code):
    pdir = parent_dir(code)
    return '%(pdir)s/samples/true-pdf.npy' % dict(pdir=pdir)

def sample_fname(code, sample_size, i):
    pdir = parent_dir(code)
    return '%(pdir)s/samples/%(num)05d/data-%(num)05d-%(i)03d.csv' % dict(pdir=pdir, num=sample_size, i=i)

def bandwidth_fname(code, sample_size):
    pdir = parent_dir(code)
    ensure_dir('%(pdir)s/kde_bw/%(num)05d/' % dict(pdir=pdir, num=sample_size))
    return '%(pdir)s/kde_bw/%(num)05d/bw.npy' % dict(pdir=pdir, num=sample_size)

# STEPS/<dist>/kde_ise_hd/<num>/ise-<start>-<block_size>.csv
def ise_hd_fname(dist_code, sample_size, start, block_size, symmetric):
    pdir = parent_dir(dist_code)
    kde_data = 'kde_sym' if symmetric else 'kde_ise_hd'
    ensure_dir('%(pdir)s/%(kde_data)s/%(num)05d/' % dict(pdir=pdir, num=sample_size, kde_data=kde_data))
    return '%(pdir)s/%(kde_data)s/%(num)05d/kde-%(block_size)03d-%(start)03d.csv' % dict(
        pdir=pdir, num=sample_size, start=start, block_size=block_size, kde_data=kde_data)

def read_true_pdf(code):
    return np.load(true_pdf_fname(code))




