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
from steps.common import sample_fname, bandwidth_fname, ise_hd_fname, read_true_pdf
from scripts2d.utils import mise_mesh

def calc_kde_ise_hd(kde, true_pdf):
    def l2_norm(v1_lin, v2_grid):
        v1_lin = v1_lin.reshape(v2_grid.shape)
        diff = v1_lin - v2_grid
        # Ok, because we know domain [0,1]x[0,1] => area = 1 x 1 = 1
        err = (diff * diff).sum()
        # extreme values are zero, which do not contribute to integral, hence correction
        # in size "_ - 1".
        return err/nns

    nns = reduce(lambda x, y: (x - 1) * (y - 1), true_pdf.shape)
    d = len(true_pdf.shape)
    # reflection to emulate [0.0:1.0:(num)j, ... as many dimensions as required, then
    # generate a nxd matrix which can be used in kde.pdf(_)
    points = np.mgrid.__getitem__(tuple([slice(0.0, 1.0, num * 1j) for num in true_pdf.shape])).reshape(d, -1).T
    pred_xy = kde.pdf(points)
    ise = l2_norm(pred_xy, true_pdf)
    hd = l2_norm(np.sqrt(pred_xy), np.sqrt(true_pdf))
    return ise, hd


class KDEIseWriter(object):
    def __init__(self, dist_code, sample_size, start, block_size):
        self.fname = ise_hd_fname(dist_code, sample_size, start, block_size)
        self.fh = None
        self.t0 = datetime.now()
        print('KDEIseWriter: Generating %s' % self.fname)

    def __enter__(self):
        self.fh = open(self.fname, 'w')
        headers = 'fname, dist_code, wave_code, n, j0, j1, k, ise, hd, elapsed_time\n'
        self.fh.write(headers)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.fh is not None:
            elapsed_time = (datetime.now() - self.t0).total_seconds()
            print('KDEIseWriter: Finished in %f secs' % elapsed_time)
            self.fh.close()

    def write_ise(self, fname, dist_code, n, ise, hd, elapsed_time):
        wave_code = 'kde'
        j0 = j1 = k = 0
        new_entry = '"%s", "%s", "%s", %d, %d, %d, %d, %f, %f, %f\n' % (
        fname, dist_code, wave_code, n, j0, j1, k, ise, hd, elapsed_time)
        ## print(new_entry)
        self.fh.write(new_entry)


def main(dist_code, sample_size, start, block_size):
    """
    calculate ISE and HD for certain samples
    :param dist_code: distribution code from utils
    :param sample_size: one of the sample sizes available
    :param start: first sample to process (starts at 0)
    :param block_size: num samples to do
    :return: None
    """
    bw_fname = bandwidth_fname(dist_code, sample_size)
    bw = np.load(bw_fname, allow_pickle=False)
    print(bw)
    ## read pdf vals for dist_code
    true_pdf = read_true_pdf(dist_code)
    with KDEIseWriter(dist_code, sample_size, start, block_size) as writer:
        for i in range(start, start + block_size):
            fname = sample_fname(dist_code, sample_size, i)
            sample = np.genfromtxt(fname, delimiter=',')
            t0 = datetime.now()
            kde = KDEMultivariate(sample, 'c' * sample.shape[1], bw=bw)
            elapsed_time = (datetime.now() - t0).total_seconds()
            ise, hd = calc_kde_ise_hd(kde, true_pdf)
            writer.write_ise(fname, dist_code, sample_size, ise, hd, elapsed_time)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step to calculate KDE for samples and corresponding ISE and HD")
    parser.add_argument('code', help='code of distribution')
    parser.add_argument('n', help='sample size', type=int)
    parser.add_argument('start', help='Start sample (0 based)', type=int)
    parser.add_argument('block_size', help='Num samples to do', type=int)
    args = parser.parse_args()

    main(args.code, args.n, args.start, args.block_size)
