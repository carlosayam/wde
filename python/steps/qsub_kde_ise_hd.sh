#!/bin/bash

DIST=$1
NUM=$2

declare -a START
if [ $DIST == mul3 ]; then
    if [ $NUM -le 512 ]; then
      START=`seq 0 128 499`
      BLOCK=128
    else
      START=`seq 0 32 499`
      BLOCK=32
    fi
else
    if [ $NUM -le 512 ]; then
      START[0]=0
      BLOCK=512
    else
      START=`seq 0 128 499`
      BLOCK=128
    fi
fi

for start in ${START[*]}; do
    cat << EOF | qsub -
#PBS -N KDE_${DIST}_${NUM}_$start
#PBS -l nodes=1:ppn=1
#PBS -l vmem=4gb
#PBS -l walltime=2:30:00
#PBS -j oe
#PBS -M z3403159@student.unsw.edu.au
#PBS -m ae

module purge
module add python/3.5.2

SW_DIR="\$PBS_O_HOME/WDE/wde/python"
. \$SW_DIR/wdeenv3/bin/activate
cd \$SW_DIR
export PYTHONPATH=.
python steps/kde_ise_hd.py $DIST $NUM $start $BLOCK
EOF
    echo "SUBMITTED steps/kde_ise_hd.py $DIST $NUM $start $BLOCK"
done
