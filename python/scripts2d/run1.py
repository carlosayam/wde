import argparse
import random
import datetime
import scripts2d.utils as u
import numpy as np
from wde.estimator import WaveletDensityEstimator
from wde.simple_estimator import SimpleWaveletDensityEstimator
from wde.estimator_via_likelihood import WaveletDensityEstimatorByLikelihood
from steps.common import grid_points, calc_true_pdf

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm


def generate_plot(dist):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.view_init(elev=0, azim=-45)
    ax.set_zlim(-2,20)
    X = np.linspace(0.0,1.0, num=256)
    Y = np.linspace(0.0,1.0, num=256)
    XX, YY = np.meshgrid(X, Y)
    Z = dist.pdf((XX, YY))
    # see http://mpastell.com/2013/05/02/matplotlib_colormaps/
    surf = ax.plot_surface(XX, YY, Z, edgecolors='k', linewidth=0.5, cmap=cm.get_cmap('BuGn'))
    #ax.set_zlim(0, 5)
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    plt.show()
    plt.close(fig)

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
    if args.true:
        generate_plot(dist)
        return
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

    ise, hd = calc_mle_ise_hd(wde.pdf.dim, wde, true_pdf)
    print(args)
    print('time', elapsed_time)
    print('ISE =', ise, 'HD = ', hd)
    generate_plot(wde)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="single wde estimator")
    parser.add_argument('code', help='code of distribution')
    parser.add_argument('wave', help='wavelet basis')
    parser.add_argument('-j0', help='initial resolution', type=int, default=0)
    parser.add_argument('-j1', help='finale resolution', type=int)
    parser.add_argument('-k', help='k in k-NN', type=int, default=1)
    parser.add_argument('-n', help='number of samples', type=int, default=1000)
    parser.add_argument('--seed', help='seed', type=int, default=random.randint(0, 999999))
    parser.add_argument('--true', help='Just plot true density', action='store_true')
    args = parser.parse_args()

    main(args)
