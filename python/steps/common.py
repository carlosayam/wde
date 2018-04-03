import os
import numpy as np
import math
import sqlite3
from scripts2d.utils import dist_from_code

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

def dbname(dist_name):
    path_name = 'data/STEPS/%s/db' % dist_name
    ensure_dir(path_name)
    return '%s/data.db' % path_name

def read_true_pdf(code):
    return np.load(true_pdf_fname(code))


def grid_points(dim, grid_n=None):
    if grid_n is None:
        grid_n = 256 if dim == 2 else 32
    points = np.mgrid.__getitem__(tuple([slice(0.0, 1.0,  grid_n * 1j) for num in range(dim)])).reshape(dim, -1).T
    return points

def calc_true_pdf(dist_code, points=None):
    dist = dist_from_code(dist_code)
    points = grid_points(dist.dim, points)
    return dist.dim, dist.pdf(points)


def create_table(conn):
    # algorithm SPW, CLW, KDE
    # wave_code, db10 (for SPW, CLW) or ''
    sql = """
    CREATE TABLE IF NOT EXISTS results (
     fname varchar(256) NOT NULL,
     algorithm varchar(8) NOT NULL,
     wave_code varchar(8) NOT NULL,
     n integer NOT NULL,
     j0 integer NOT NULL,
     j1 integer NOT NULL,
     k integer NOT NULL,
     ise real NOT NULL,
     hd real NOT NULL,
     etime real NOT NULL
     )
    """
    conn.execute(sql)
    print('results created')

def connect(dist_name, create=False):
    fname_db = dbname(dist_name)
    if create and not os.path.isfile(fname_db):
        print('Creating %s' % fname_db)
        conn = sqlite3.connect(fname_db)
        create_table(conn)
    else:
        print('Opening %s' % fname_db)
        conn = sqlite3.connect(fname_db)
    return conn

def l2_norm(v1_lin, v2_lin, nns):
    diff = v1_lin - v2_lin
    # Ok, because we know domain [0,1]x[0,1] => area = 1 x 1 = 1
    err = (diff * diff).sum()
    # extreme values are zero, which do not contribute to integral, hence correction
    # in size "_ - 1".
    return err/nns

def calc_spwe_ise(dim, wde, true_pdf, points=None):
    points = grid_points(dim, points)
    pred_xy = wde.pdf(points)
    nns = points.shape[0]
    ise = l2_norm(pred_xy, true_pdf, nns)
    return ise

def calc_spwe_ise_hd(dim, wde, true_pdf):
    points = grid_points(dim)
    pred_xy = wde.pdf(points)
    nns = points.shape[0]
    ise = l2_norm(pred_xy, true_pdf, nns)
    hd = l2_norm(np.sqrt(pred_xy), np.sqrt(true_pdf), nns)
    return ise, hd