import argparse
import random
import datetime
import scripts2d.utils as u
import numpy as np


def main(dist_code):
    # for each dir_n in data/RESP/dist_code/data2d/samples/*
    #   mkdir dir_n/kde
    #   for each fname in dir_n/*.csv
    #     read sample fname
    #     calculate best kde w/ cross validation
    #     write bandwidth, time to dir_n/kde_bw/fname.bw
    #   cat kde_bw/
    #   launch calculate_ise for dir_n
    #   launch calculate_hd for dir_n

    # - calculate_ise (dir_n)
    # for each sample in dir_n:
    #    read sample
    #    read bandwidth
    #    create estimator
    #    calculate ise
    #    write to kde_ise/fname.tmp
    # cat kde_ise/*.tmp > kde_ise/all.csv
    # remove kde_ise/*.tmp

    # - calculate_hd (dir_n)
    # for each sample in dir_n:
    #    read sample
    #    read bandwidth
    #    create estimator
    #    calculate hd
    #    write to kde_hd/fname.tmp
    # cat kde_hd/*.tmp > kde_hd/all.csv
    # remove kde_hd/*.tmp

    dist = u.dist_from_code(args.code)
    wave_code = args.wave
    random.seed(args.seed)
    data = dist.rvs(args.n)
    true_z = dist.pdf(u.mise_mesh(dist.dim))
    t0 = datetime.datetime.now()
    wde.fit(data)
    elapsed_time = (datetime.datetime.now() - t0).total_seconds()
    pred_z = wde.pdf(u.mise_mesh(d=dist.dim))
    ise = u.calc_ise(wde.pdf, true_z)
    print(args)
    print('time', elapsed_time)
    print('true shape', true_z.shape, 'pred shape', pred_z.shape)
    print('calc_ise', ise)
    # KDE
    kde = stats.gaussian_kde(data.T)
    grid = u.mise_mesh(dist.dim)
    pos = np.stack(grid, axis=0)
    vals = np.moveaxis(pos, 0, -1)
    vals = vals.reshape(-1, vals.shape[-1])
    print(vals.shape)
    pred_z = kde(vals.T)
    print('grid shape', true_z.shape)
    print('flat shape', true_z.flatten().shape)
    print('kde_ise', mean_squared_error(true_z.flatten(), pred_z))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="single wde estimator")
    parser.add_argument('dist_code', help='code of distribution')
    args = parser.parse_args()

    main(args.dist_code)
