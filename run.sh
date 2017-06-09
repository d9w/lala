#!/bin/sh

HOME=/users/p16043/wilson/lala

WORK_DIR=/tmpdir/$LOGNAME/dennis/$SLURM_JOB_ID/$SLURM_TASK_PID
mkdir -p $WORK_DIR

cd $WORK_DIR
cp $HOME/cond/single.jl $WORK_DIR
#cp $HOME/res/defaults.json $WORK_DIR
#cp $HOME/res/worm_defaults.json $WORK_DIR/shape_defaults.json
#cp $HOME/res/ranges.json $WORK_DIR
#cp $HOME/julia/runcmaes.jl $WORK_DIR
#cp $HOME/julia/cmaes.jl $WORK_DIR
#cp $HOME/res/worm.shape $WORK_DIR
#cp $HOME/bin/console $WORK_DIR
#cd res
export JULIA_PKGDIR=/users/p16043/wilson/.julia
julia single.jl 0 $SLURM_TASK_PID
