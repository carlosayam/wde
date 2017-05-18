#!/usr/bin/env python

from __future__ import division

import os, sys
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

PBS="""
#!/bin/bash

#PBS -N ARRAY
#PBS -l nodes=1:ppn=1
#PBS -l vmem=1gb
#PBS -l walltime=2:00:00
#PBS -j oe
#PBS -t 1-%d

cd $PBS_O_WORKDIR
 
module purge
module add python/2.7.12

RESP_DIR="$PBS_O_WORKDIR"

mkdir -p $RESP_DIR/data2d/samples
mkdir -p $RESP_DIR/data2d/resp

cd $RESP_DIR

. wdeenv/bin/activate
./scripts2d/run_for.py 1000 $PBS_ARRAYID

"""[1:-1]

def write_pbs(num):
    fname = 'data2d/run.pbs'
    with open(fname, 'wa') as f:
        f.write(PBS % (num // 1000))
    sys.exec('chmod a+x %s' % fname) #??

def main():
    plans = []
    beta = u.dist_from_code('beta')
    multi = u.dist_from_code('mult')
    for dist in [beta, multi]:
        for ix in range(9):
            n = 1000 + ix * 2500
            for i in range(5):
                data = dist.rvs(n)
                fname = write_sample(dist, n, i, data)
                for j0 in range(0, 4):
                    for j1 in range(j0 - 1, 8):
                        k = 1
                        while k < int(sqrt(n)/4):
                            plans.append(dict(fname=fname, wave_code=dist.code, j0=j0, j1=j1, k=k))
                            k = 2 * k
    u.write_plans(plans)
    write_pbs(len(plans))
    u.write_dist(beta)
    u.write_dist(multi)
    print 'Done %d' % len(plans)

if __name__ == "__main__":
    main()
