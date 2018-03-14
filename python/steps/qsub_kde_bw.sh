#!/bin/bash
DIST=$1
NUM=$2

if [ $NUM -le 1024 ]; then
  WALL_TIME='3:30:00'
else
  WALL_TIME='8:30:00'
fi

cat << EOF | qsub -
#PBS -N KDE_BW_${DIST}_${NUM}
#PBS -l nodes=1:ppn=1
#PBS -l vmem=4gb
#PBS -l walltime=${WALL_TIME}
#PBS -j oe
#PBS -M z3403159@student.unsw.edu.au
#PBS -m ae

module purge
module add python/3.5.2

SW_DIR="\$PBS_O_HOME/WDE/wde/python"
. \$SW_DIR/wdeenv3/bin/activate
cd \$SW_DIR
export PYTHONPATH=.
python steps/kde_bandwidth.py $DIST $NUM
EOF

echo "SUBMITTED steps/kde_bandwidth.py $DIST $NUM"

