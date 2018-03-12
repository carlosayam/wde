#
# $ python steps/kde_ise_hd.py <dist> <num> <start> <block_size>
#
# read bandwidth from STEPS/<dist>/kde_bw/<num>/bw.npy as numpy array
# read sample files from STEPS/<dist>/samples/<num>/data-<num>-{000_NNN}.csv
# It reads samples starting from <start> (0 to N) and processes only
# <block_size> samples.
# It leaves the result in STEPS/<dist>/kde_ise_hd/<num>/ise-<start>-<block_size>.csv
#

import argparse
import numpy as np
import math
from datetime import datetime
from steps.common import sample_fname, parent_dir, ensure_dir, Adder
from scripts2d.utils import mise_mesh, dist_from_code
from wde.estimator import WaveletDensityEstimator

def ensure_samples_dir(dist_code, sample_size):
    pdir = parent_dir(dist_code)
    the_dir = '%(pdir)s/samples/%(num)05d' % dict(pdir=pdir, num=sample_size)
    ensure_dir(the_dir, False) ## throw error if exists
    return the_dir


def sample_fname(the_dir, sample_size, index):
    return '%(the_dir)s/data-%(num)05d-%(index)03d.csv' % \
           dict(the_dir=the_dir, num=sample_size, index=index)


def main(dist_code, sample_size):
    """
    Generates 500 samples of given size
    :param dist_code:
    :param sample_size:
    :return:
    """
    ## read pdf vals for dist_code
    dist = dist_from_code(dist_code)
    the_dir = ensure_samples_dir(dist_code, sample_size)
    t0 = datetime.now()
    for index in range(500):
        data = dist.rvs(sample_size)
        fname = sample_fname(the_dir, sample_size, index)
        np.savetxt(fname, data, fmt='%f', delimiter=',')
    elapsed_time = (datetime.now() - t0).total_seconds()
    print('Done gen_samples %s, %s in %f' % (dist_code, sample_size, elapsed_time))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step to generate samples of given size")
    parser.add_argument('code', help='code of distribution')
    parser.add_argument('n', help='sample size', type=int)
    args = parser.parse_args()

    main(args.code, args.n)
