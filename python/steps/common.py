import os
import numpy as np

if 'PBS_O_HOME' in os.environ:
    PBS_O_HOME = os.environ['PBS_O_HOME']
else:
    PBS_O_HOME = 'data'


def ensure_dir(directory):
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except FileExistsError as ex:
            # ignore race condition between multiple jobs
            pass

def parent_dir(code):
    return '%(home)s/STEPS/%(dist_code)s/' % dict(home=PBS_O_HOME, dist_code=code)

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
def ise_hd_fname(dist_code, sample_size, start, block_size):
    pdir = parent_dir(dist_code)
    ensure_dir('%(pdir)s/kde_ise_hd/%(num)05d/' % dict(pdir=pdir, num=sample_size))
    return '%(pdir)s/kde_ise_hd/%(num)05d/kde-%(block_size)03d-%(start)03d.csv' % dict(
        pdir=pdir, num=sample_size, start=start, block_size=block_size)

def read_true_pdf(code):
    return np.load(true_pdf_fname(code))




