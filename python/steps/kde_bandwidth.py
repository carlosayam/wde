#
# $ python steps/kde_bandwidth.py <dist> <num>
#
# from STEPS/<dist>/<num>/samples pick first 10 samples
# calculate bandwitdh
# store bandwidth into STEPS/<dist>/kde_bw/<num>/bw.npy as numpy array
#

import argparse
import numpy as np
from statsmodels.nonparametric.kernel_density import KDEMultivariate
from steps.common import sample_fname, bandwidth_fname

def main(code, sample_size):
    ## calculate average bw for <num> samples
    bw = None
    num = 10
    for i in range(num):
        fname = sample_fname(code, sample_size, i)
        sample = np.genfromtxt(fname, delimiter=',')
        kde = KDEMultivariate(sample, 'c' * sample.shape[1], bw='cv_ml')
        if bw is None:
            bw = kde.bw
        else:
            bw += kde.bw
    bw_fname = bandwidth_fname(code, sample_size)
    bw = bw / num # average
    np.save(bw_fname, bw, allow_pickle=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="kde bandwidth calculator")
    parser.add_argument('code', help='code of distribution')
    parser.add_argument('n', help='sample size', type=int)
    args = parser.parse_args()

    main(args.code, args.n)
