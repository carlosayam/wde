#
# $ python steps/kde_bandwidth.py <dist> <num>
#
# from STEPS/<dist>/<num>/samples pick first 10 samples
# calculate bandwitdh
# store bandwidth into STEPS/<dist>/kde_bw/<num>/bw.npy as numpy array
#
import sys
import argparse
from datetime import datetime
import numpy as np
from wde.estimator import WaveletDensityEstimator
from steps.common import sample_fname, parent_dir, ensure_dir, calc_spwe_ise, calc_true_pdf
from wde.thresholding import soft_threshold

POINTS=64

def threshold_fname(code, wave_code, kind):
    pdir = parent_dir(code)
    ensure_dir('%(pdir)s/spwe_th' % dict(pdir=pdir))
    return '%(pdir)s/spwe_th/%(wave_code)_%(kind)s.npy' % dict(pdir=pdir, wave_code=wave_code, kind=kind)

def soft_threshold_calc(wde, true_pdf):
    # soft thresholding calculation
    tt0 = 0
    tt1 = 1
    err = 1
    best_ise = calc_spwe_ise(wde.dim, wde, true_pdf, POINTS)
    best_tt = tt1
    best_i = 7
    while err > 0.00001:
        print('.', end='')
        sys.stdout.flush()
        err = 0
        for i in range(8):
            tt_i = tt0 + (tt1 - tt0) * i / 7
            wde.thresholding = soft_threshold(tt_i)
            wde.pdf = wde.calc_pdf()
            new_ise = calc_spwe_ise(wde.dim, wde, true_pdf, POINTS)
            if new_ise < best_ise:
                err = best_ise - new_ise
                best_ise = new_ise
                best_tt = tt_i
                best_i = i
        i0 = max(best_i - 2, 0)
        i1 = min(best_i + 2, 7)
        tt0, tt1 = tt0 + (tt1 - tt0) * i0 / 7, tt0 + (tt1 - tt0) * i1 / 7
    print('.')
    return best_tt

def main(code, wave_code, sample_size, j0, j1):
    ## calculate average threshold for based on (first) 50 samples
    bw = None
    num = 50
    t0 = datetime.now()
    dim, true_pdf = calc_true_pdf(code, POINTS)
    thresholds = []
    t0 = datetime.now()
    for i in range(num):
        fname = sample_fname(code, sample_size, i)
        print('+', end='')
        sample = np.genfromtxt(fname, delimiter=',')
        wde = WaveletDensityEstimator(wave_code, k=1, j0=j0, j1=j1)
        wde.fit(sample)
        thresholds.append(soft_threshold_calc(wde, true_pdf))
    fname = threshold_fname(code, wave_code, 'soft')
    threshold = np.array(thresholds).mean()
    np.save(fname, threshold, allow_pickle=False)
    elapsed_time = (datetime.now() - t0).total_seconds()
    print('\nDone soft threshold %s, %s in %f' % (code, sample_size, elapsed_time))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="spwe threshold calculator")
    parser.add_argument('code', help='code of distribution')
    parser.add_argument('wave', help='wavelt code to use')
    parser.add_argument('n', help='sample size', type=int, default=1024)
    parser.add_argument('j0', help='initial level', type=int, default=3)
    parser.add_argument('j1', help='final level', type=int, default=6)
    args = parser.parse_args()

    main(args.code, args.wave, args.n, args.j0, args.j1)
