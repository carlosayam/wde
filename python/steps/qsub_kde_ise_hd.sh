#!/bin/bash

$DIST=$1
$NUM=$2

declare -a START
if [ $NUM <= 512 ]; then
  START[0]=0
  BLOCK=512
else
  START=(0, 128, 256, 384)
  BLOCK=128
fi

for start in ${START[*]}; do
    cat << EOF > tmp.pbs
    #PBS -N PLAN_$DIST_$NUM_$start
    #PBS -l nodes=1:ppn=1
    #PBS -l vmem=4gb
    #PBS -l walltime=1:00:00
    #PBS -j oe
    #PBS -M z3403159@student.unsw.edu.au
    #PBS -m ae

    module purge
    module add python/3.5.2

    SW_DIR="$PBS_O_HOME/WDE/wde/python"

    $SW_DIR/steps/kde_ise_hd.py $DIST $NUM $start $BLOCK
EOF

    qsub tmp.pbs
    echo "SUBMITTED steps/kde_ise_hd.py $DIST $NUM $start $BLOCK"
    sleep 10
done
