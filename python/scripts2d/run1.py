import argparse
import random
import datetime
import scripts2d.utils as u
import numpy as np
from wde.estimator import WaveletDensityEstimator
from wde.simple_estimator import SimpleWaveletDensityEstimator
from wde.estimator_via_likelihood import WaveletDensityEstimatorByLikelihood
from sklearn.metrics import mean_squared_error
from scipy import stats
from steps.common import grid_points, calc_true_pdf

def calc_mle_ise_hd(dim, wmle, true_pdf):
    def l2_norm(v1_lin, v2_lin):
        diff = v1_lin - v2_lin
        # Ok, because we know domain [0,1]x[0,1] => area = 1 x 1 = 1
        err = (diff * diff).sum()
        # extreme values are zero, which do not contribute to integral, hence correction
        # in size "_ - 1".
        return err/nns

    points = grid_points(dim)
    pred_xy = wmle.pdf(points)
    nns = points.shape[0]
    ise = l2_norm(pred_xy, true_pdf)
    hd = l2_norm(np.sqrt(pred_xy), np.sqrt(true_pdf))
    return ise, hd


def main(args):
    dist = u.dist_from_code(args.code)
    wave_code = args.wave
    random.seed(args.seed)
    data = dist.rvs(args.n)
    true_z = dist.pdf(u.mise_mesh(dist.dim))
    if wave_code[0:4] == 'sim-':
        wde = SimpleWaveletDensityEstimator(wave_code[4:], j0 = args.j0, j1 = args.j1)
    elif wave_code[0:4] == 'mle-':
        wde = WaveletDensityEstimatorByLikelihood(wave_code[4:], j0=args.j0, j1=args.j1)
    else:
        wde = WaveletDensityEstimator(wave_code, k = args.k, j0 = args.j0, j1 = args.j1)
    t0 = datetime.datetime.now()
    wde.fit(data)
    elapsed_time = (datetime.datetime.now() - t0).total_seconds()
    dim, true_pdf = calc_true_pdf(args.code)

    pred_z = wde.pdf(u.mise_mesh(d=dist.dim))
    ise, hd = calc_mle_ise_hd(wde.pdf.dim, wde, true_pdf)
    print(args)
    print('time', elapsed_time)
    print('ISE =', ise, 'HD = ', hd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="single wde estimator")
    parser.add_argument('code', help='code of distribution')
    parser.add_argument('wave', help='wavelet basis')
    parser.add_argument('-j0', help='initial resolution', type=int, default=0)
    parser.add_argument('-j1', help='finale resolution', type=int)
    parser.add_argument('-k', help='k in k-NN', type=int, default=1)
    parser.add_argument('-n', help='number of samples', type=int, default=1000)
    parser.add_argument('--seed', help='seed', type=int, default=random.randint(0, 999999))
    args = parser.parse_args()

    main(args)
