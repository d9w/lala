#!/bin/sh

# requires env variable $JOB_ID

HOME=/users/p16043/wilson/lala
WORK_DIR=/tmpdir/wilson/dennis/$JOB_ID/
RESULTS_DIR=$HOME/results/$JOB_ID/

mkdir -p $RESULTS_DIR

i=0
for proc in $WORK_DIR/*
do
    # cat $proc/log | grep 'R:' | rev | cut -d ':' -f 1 | rev > $RESULTS_DIR/train_$i.log
    cp $proc/cmaes.log $RESULTS_DIR/$i.log
    echo $i
    let i=i+1
done
