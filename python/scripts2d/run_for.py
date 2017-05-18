#!/usr/bin/env python

import os, sys
import pandas

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')

from wde.estimator import WaveletDensityEstimator
import scripts2d.utils as u

# parameters
# dist : B = beta, M = multi
# bag_size : 1000
# bag_number
# - reads plan.csv and runs all data from bag_size * (bag_number - 1) .. + bag_size
# -- produces files in
# data2d/results/{sample_file}-pdf.csv : generated PDF in grid/numpy shape
# data2d/results/{sample_file}-mise.csv :  calculated MISE

def exec_plan(row):
    fname, wave_name, j0, j1, k = row['fname'], row['wave_code'], row['j0'], row['j1'], row['k']
    data, dist = u.read_sample(fname)
    wde = WaveletDensityEstimator(wave_name, k = k, j0 = j0, j1 = j1)
    wde.fit(data)
    mise = u.calc_mise(wde.pdf, dist.pdf)

