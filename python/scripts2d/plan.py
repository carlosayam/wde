#!/usr/bin/env python
from __future__ import division

import math
import os
import random
import sys
import numpy as np
import pandas

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')

import scripts2d.utils as u

# produces files in
# data2d/samples/wde.{dist}.{n}-{i}.csv : generated samples
# data2d/plan.csv : generated plan
# data2d/run_all.pbs : PBS plan

# dist = B (beta) | M (multinomial)
# n = 1000..21000, step=2500
# i = 500 samples of this size

# for each sample, (4000 per dist)
#   gen item plan for j0=0..2 (incl), j1=j0-..5, k=1,2,4,8,..sqrt(n)/4
# generate then PBS that runs "run_for" in parallel #plans nodes

BAG_SIZE=1000

PBS="""
#!/bin/bash

#PBS -N ARRAY
#PBS -l nodes=1:ppn=1
#PBS -l vmem=1gb
#PBS -l walltime=2:00:00
#PBS -j oe
#PBS -t 1-\%d

module purge
module add python/2.7.12

RESP_DIR="$PBS_O_INITDIR/RESP"
SW_DIR="$PBS_O_INITDIR/WDE/wde/python"

mkdir -p $RESP_DIR
cd $RESP_DIR

cd $RESP_DIR

. $SW_DIR/wdeenv/bin/activate
$SW_DIR/scripts2d/run_for.py %d $PBS_ARRAYID
"""[1:] % BAG_SIZE

def write_pbs(num):
    fname = 'data2d/run.pbs'
    with open(fname, 'wa') as f:
        f.write(PBS % (num // BAG_SIZE + 1))

def gen_samples():
    beta = u.dist_from_code('beta')
    multi = u.dist_from_code('mult')
    for dist in [beta, multi]:
        for ix in range(3): # 9
            n = 100 + ix * 250 # 1000, 2500
            for i in range(2): # 500
                data = dist.rvs(n)
                fname = u.write_sample(dist, n, i, data)
                yield fname, n, dist

def main():
    wave_code = sys.argv[3]
    plans = []
    for fname, n, dist in gen_samples():
        for wave_code in ['db1', 'db3', 'db5', 'sym4', 'coif1']:
            for j0 in range(0, 4):
                for j1 in range(j0, 8):
                    k = 1
                    while k < int(math.sqrt(n)/4):
                        plans.append(dict(fname=fname, dist_code=dist.code, wave_code=wave_code, j0=j0, j1=j1, k=k, rand=random.random()))
                        k = 2 * k
    u.write_plans(plans)
    u.write_dist_pdf(u.dist_from_code('beta'))
    u.write_dist_pdf(u.dist_from_code('mult'))
    u.empty_ise()
    write_pbs(len(plans))
    print 'Done %d' % len(plans)

if __name__ == "__main__":
    main()
