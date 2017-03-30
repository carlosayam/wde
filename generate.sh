#!/bin/bash

for j1 in 4 5; do
    let nsamples=$((j1 * 1000 + 1000))
    cat generate.pbs | sed "s/_J1_/$j1/" | sed "s/_NSAMPLES_/$nsamples/" > temp.pbs
    qsub temp.pbs
    sleep 10
done
