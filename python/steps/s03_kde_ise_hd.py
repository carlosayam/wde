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
from functools import reduce
from datetime import datetime
from statsmodels.nonparametric.kernel_density import KDEMultivariate
from steps.common import sample_fname, bandwidth_fname, ise_hd_fname, Adder, grid_points, calc_true_pdf


def calc_kde_ise_hd(dim, kde, true_pdf):
    def l2_norm(v1_lin, v2_lin):
        diff = v1_lin - v2_lin
        # Ok, because we know domain [0,1]x[0,1] => area = 1 x 1 = 1
        err = (diff * diff).sum()
        # extreme values are zero, which do not contribute to integral, hence correction
        # in size "_ - 1".
        return err/nns

    points = grid_points(dim)
    pred_xy = kde.pdf(points)
    nns = points.shape[0]
    ise = l2_norm(pred_xy, true_pdf)
    hd = l2_norm(np.sqrt(pred_xy), np.sqrt(true_pdf))
    return ise, hd


class KDEIseWriter(object):
    def __init__(self, dist_code, sample_size, start, block_size, symmetric, dry=False):
        self.fname = ise_hd_fname(dist_code, sample_size, start, block_size, symmetric)
        self.fh = None
        self.t0 = datetime.now()
        self.num = block_size
        self.start = start
        # quick average
        self.dry = dry
        self.ise_adder = Adder('ISE')
        self.hd_adder = Adder('HD')
        self.etime_adder = Adder('E.TIME')
        print('KDEIseWriter: Generating %s' % self.fname)

    def __enter__(self):
        self.fh = open(self.fname, 'w', 1) if not self.dry else None
        if self.start == 0 and not self.dry:
            headers = 'fname, dist_code, wave_code, n, j0, j1, k, ise, hd, elapsed_time\n'
            self.fh.write(headers)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        elapsed_time = (datetime.now() - self.t0).total_seconds()
        print('KDEIseWriter: Finished in %f secs' % elapsed_time)
        print('Stats: ', self.ise_adder, self.hd_adder, self.etime_adder)
        if self.fh is not None:
            self.fh.close()
        if exc_type is not None:
            if self.num <= 0:
                ## block exception, we have done all
                return True

    def write_ise(self, fname, dist_code, n, ise, hd, elapsed_time):
        wave_code = 'kde'
        j0 = j1 = k = 0
        new_entry = '"%s", "%s", "%s", %d, %d, %d, %d, %f, %f, %f\n' % (
        fname, dist_code, wave_code, n, j0, j1, k, ise, hd, elapsed_time)
        #print(new_entry)
        if not self.dry:
            self.fh.write(new_entry)
        self.ise_adder.add_v(ise)
        self.hd_adder.add_v(hd)
        self.etime_adder.add_v(elapsed_time)
        self.num -= 1


def main(dist_code, sample_size, start, block_size, symmetric):
    """
    calculate ISE and HD for certain samples
    :param dist_code: distribution code from utils
    :param sample_size: one of the sample sizes available
    :param start: first sample to process (starts at 0)
    :param block_size: num samples to do
    :param symmetric: use same bandwidth in all dimensions
    :return: None
    """
    bw_fname = bandwidth_fname(dist_code, sample_size)
    bw = np.load(bw_fname, allow_pickle=False)
    if symmetric:
        bw = np.ones(bw.shape[0]) * np.mean(bw)
    print('Bandwidth', bw, '- Symmetric', symmetric)
    ## read pdf vals for dist_code
    dim, true_pdf = calc_true_pdf(dist_code)
    with KDEIseWriter(dist_code, sample_size, start, block_size, symmetric) as writer:
        for i in range(start, start + block_size):
            if i >= 500:
                continue
            fname = sample_fname(dist_code, sample_size, i)
            sample = np.genfromtxt(fname, delimiter=',')
            t0 = datetime.now()
            kde = KDEMultivariate(sample, 'c' * sample.shape[1], bw=bw)
            elapsed_time = (datetime.now() - t0).total_seconds()
            ise, hd = calc_kde_ise_hd(dim, kde, true_pdf)
            writer.write_ise(fname, dist_code, sample_size, ise, hd, elapsed_time)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step to calculate KDE for samples and corresponding ISE and HD")
    parser.add_argument('code', help='code of distribution')
    parser.add_argument('n', help='sample size', type=int)
    parser.add_argument('start', help='Start sample (0 based)', type=int)
    parser.add_argument('block_size', help='Num samples to do', type=int)
    parser.add_argument('--symmetric', help='same bandwidth all directions', action='store_true')
    args = parser.parse_args()

    main(args.code, args.n, args.start, args.block_size, args.symmetric)
