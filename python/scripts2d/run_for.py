#!/usr/bin/env python
from __future__ import division

import datetime
import os
import sys
import pandas

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')

from wde.estimator import WaveletDensityEstimator
import scripts2d.utils as u

# parameters
# bag_size : 1000
# bag_number
# - reads plan.csv and runs all data from bag_size * (bag_number - 1) .. + bag_size
# -- produces files in
# data2d/results/{sample_file}-pdf.csv : generated PDF in grid/numpy shape
# data2d/results/{sample_file}-mise.csv :  calculated MISE

def exec_plan(fh_ise, fh_coeffs, row):
    fname, wave_code, dist_code, j0, j1, k, do_coeffs = row['fname'], row['wave_code'], row['dist_code'], row['j0'], row['j1'], row['k'], row['do_coeffs']
    data = u.read_sample(fname)
    n = len(data)
    pdf_vals = u.read_dist_pdf()
    wde = WaveletDensityEstimator(wave_code, k = k, j0 = j0, j1 = j1)
    t0 = datetime.datetime.now()
    wde.fit(data)
    elapsed_time = (datetime.datetime.now() - t0).total_seconds()
    ise = u.calc_ise(wde.pdf, pdf_vals)
    u.write_ise(fh_ise, fname, dist_code, wave_code, n, j0, j1, k, ise, elapsed_time)
    if do_coeffs:
        u.write_coeffs(fh_coeffs, fname, dist_code, wave_code, n, j0, j1, k, wde)

def main():
    bag_size = int(sys.argv[1])
    bag_number = int(sys.argv[2])
    plans = u.read_plans(bag_size, bag_number)
    with open("data2d/ise-%04d.csv" % bag_number, "w") as fh_ise:
        with open("data2d/coeffs-%04d.csv" % bag_number, "w") as fh_coeffs:
            if bag_number == 1:
                headers = 'fname, dist_code, wave_code, n, j0, j1, k, ise, elapsed_time\n'
                fh_ise.write(headers)
                headers = 'fname, dist_code, wave_code, n, j0, j1, k, j, qx, zs, coeff\n'
                fh_coeffs.write(headers)
            for _, row in plans.iterrows():
                exec_plan(fh_ise, fh_coeffs, row)


if __name__ == "__main__":
    main()
