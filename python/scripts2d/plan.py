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

BAG_SIZE=500

PBS="""
#!/bin/bash

#PBS -N ARRAY
#PBS -l nodes=1:ppn=1
#PBS -l vmem=2gb
#PBS -l walltime=5:30:00
#PBS -j oe

#PBS -t 1-%%d

module purge
module add python/2.7.12

RESP_DIR="$PBS_O_HOME/RESP/%s/%s"
SW_DIR="$PBS_O_HOME/WDE/wde/python"

mkdir -p $RESP_DIR
cd $RESP_DIR

. $SW_DIR/wdeenv/bin/activate
$SW_DIR/scripts2d/run_for.py %d $PBS_ARRAYID
"""[1:]

def write_pbs(dist_code, wave_code, num):
    fname = 'data2d/run.pbs'
    pbs = PBS % (dist_code, wave_code, BAG_SIZE)
    with open(fname, 'w') as f:
        f.write(pbs % (num // BAG_SIZE + 1))

def ns_for(dist):
    if dist.code == 'mul2':
        return [64,96,128,192,256,384]
    return [128,256,512,1024,2048,4096]

def gen_samples(dist):
    for n in ns_for(dist):
        for i in range(500):
            data = dist.rvs(n)
            fname = u.write_sample(n, i, data)
            yield fname, n

def main(dist_code, wave_code):
    plans = []
    dist = u.dist_from_code(dist_code)
    for fname, n in gen_samples(dist):
        for j0 in range(0, 4):
            j1 = j0 - 1 # single level
            if wave_code[0:4] == 'sim-':
                plans.append(dict(fname=fname, dist_code=dist.code, wave_code=wave_code, j0=j0, j1=j1, k=1, rand=random.random()))
            else:
                k = 1
                while k * k < n:
                    plans.append(dict(fname=fname, dist_code=dist.code, wave_code=wave_code, j0=j0, j1=j1, k=k, rand=random.random()))
                    k = 2 * k
    u.write_plans(plans)
    u.write_dist_pdf(u.dist_from_code(dist_code))
    write_pbs(dist_code, wave_code, len(plans))
    print 'Done %d' % len(plans)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
