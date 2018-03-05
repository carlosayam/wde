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

def grid_points(dim):
    grid_n = 256 if dim == 2 else 64
    points = np.mgrid.__getitem__(tuple([slice(0.0, 1.0,  grid_n * 1j) for num in range(dim)])).reshape(dim, -1).T
    return points


def calc_true_pdf(dist_code):
    dist = dist_from_code(dist_code)
    grid_n = 256 if dist.dim == 2 else 64
    points = grid_points(dist.dim)
    return dist.dim, dist.pdf(points)


def calc_spwe_ise_hd(dim, wde, true_pdf):
    def l2_norm(v1_lin, v2_lin):
        diff = v1_lin - v2_lin
        # Ok, because we know domain [0,1]x[0,1] => area = 1 x 1 = 1
        err = (diff * diff).sum()
        # extreme values are zero, which do not contribute to integral, hence correction
        # in size "_ - 1".
        return err/nns

    points = grid_points(dim)
    pred_xy = wde.pdf(points)
    nns = points.shape[0]
    ise = l2_norm(pred_xy, true_pdf)
    hd = l2_norm(np.sqrt(pred_xy), np.sqrt(true_pdf))
    return ise, hd

def spwe_ise_hd_fname(dist_code, wave, j0, j1, k, sample_size, start, block_size):
    pdir = parent_dir(dist_code)
    j1str = '__' if j1 == None else ('%02d' % j1)
    the_dir = '%(pdir)s/spwe/%(wave)s/wde-%(j0)02d-%(j1str)s-%(k)03d/%(num)05d' % \
               dict(pdir=pdir, wave=wave, j0=j0, j1str=j1str, k=k, num=sample_size)
    ensure_dir(the_dir)
    return '%(the_dir)s/ise-%(block_size)03d-%(start)03d.csv' % \
           dict(the_dir=the_dir, block_size=block_size, start=start)


class SPWEIseWriter(object):
    def __init__(self, dist_code, wave, j0, j1, k, sample_size, start, block_size):
        self.fname = spwe_ise_hd_fname(dist_code, wave, j0, j1, k, sample_size, start, block_size)
        self.dist_code = dist_code
        self.wave = wave
        self.j0 = j0
        if j1 == None:
            self.j1 = j0 - 1
        else:
            self.j1 = j1
        self.k = k
        self.n = sample_size
        self.fh = None
        self.t0 = datetime.now()
        self.num = block_size
        self.start = start
        # quick average
        self.ise_adder = Adder('ISE')
        self.hd_adder = Adder('HD')
        self.etime_adder = Adder('E.TIME')
        print('SPWEIseWriter: Generating %s' % self.fname)

    def __enter__(self):
        self.fh = open(self.fname, 'w', 1)
        if self.start == 0:
            headers = 'fname, dist_code, wave_code, n, j0, j1, k, ise, hd, elapsed_time\n'
            self.fh.write(headers)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.fh is not None:
            elapsed_time = (datetime.now() - self.t0).total_seconds()
            print('SPWEIseWriter: %s done in %f secs' % (self.fname, elapsed_time))
            print('Stats: ', self.ise_adder, self.hd_adder, self.etime_adder)
            self.fh.close()
        if exc_type is not None:
            if self.num <= 0:
                ## block exception, we have done all
                return True

    def write_ise(self, fname, ise, hd, elapsed_time):
        new_entry = '"%s", "%s", "%s", %d, %d, %d, %d, %f, %f, %f\n' % \
                    (fname, self.dist_code, self.wave, self.n, self.j0, self.j1,\
                     self.k, ise, hd, elapsed_time)
        #print(new_entry)
        self.fh.write(new_entry)
        self.ise_adder.add_v(ise)
        self.hd_adder.add_v(hd)
        self.etime_adder.add_v(elapsed_time)
        self.num -= 1


def main(dist_code, wave, j0, j1, k, sample_size, start, block_size):
    """
    Calculate ISE and HD for Shape-preserving estimator for certain samples
    :param dist_code:
    :param wave:
    :param j0:
    :param j1:
    :param k:
    :param sample_size:
    :param start:
    :param block_size:
    :return:
    """
    ## read pdf vals for dist_code
    dim, true_pdf = calc_true_pdf(dist_code)
    with SPWEIseWriter(dist_code, wave, j0, j1, k, sample_size, start, block_size) as writer:
        for i in range(start, start + block_size):
            if i >= 500:
                continue
            fname = sample_fname(dist_code, sample_size, i)
            sample = np.genfromtxt(fname, delimiter=',')
            t0 = datetime.now()
            wde = WaveletDensityEstimator(wave, k = k, j0 = j0, j1 = j1)
            wde.fit(sample)
            elapsed_time = (datetime.now() - t0).total_seconds()
            ise, hd = calc_spwe_ise_hd(dim, wde, true_pdf)
            writer.write_ise(fname, ise, hd, elapsed_time)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step to calculate SP for samples and corresponding ISE and HD")
    parser.add_argument('code', help='code of distribution')
    parser.add_argument('wave', help='code of wavelet')
    parser.add_argument('j0', help='j_0', type=int)
    parser.add_argument('k', help='j_0', type=int)
    parser.add_argument('n', help='sample size', type=int)
    parser.add_argument('--start', help='Start sample (0 based)', type=int, default=0)
    parser.add_argument('--block_size', help='Num samples to do', type=int, default=500)
    parser.add_argument('--j1', help='j_1', type=int, default=None)
    args = parser.parse_args()

    main(args.code, args.wave, args.j0, args.j1, args.k, args.n, args.start, args.block_size)
