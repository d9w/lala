#!/bin/sh

HOME=/users/p16043/wilson/lala

WORK_DIR=/tmpdir/$LOGNAME/dennis/$SLURM_JOB_ID/$SLURM_TASK_PID
mkdir -p $WORK_DIR

cd $WORK_DIR
cp $HOME/config/defaults.json $WORK_DIR
cp $HOME/config/stingray_defaults.json $WORK_DIR/shape_defaults.json
cp $HOME/config/ranges.json $WORK_DIR
cp $HOME/julia/runcmaes.jl $WORK_DIR
cp $HOME/julia/cmaes.jl $WORK_DIR
cp $HOME/config/stingray.shape $WORK_DIR
cp $HOME/config/stingray_neurons.shape $WORK_DIR
cp $HOME/bin/console $WORK_DIR
export JULIA_PKGDIR=/users/p16043/wilson/.julia
julia runcmaes.jl $SLURM_TASK_PID
