#!/bin/sh

# requires env variable $JOB_ID

HOME=/users/p16043/wilson/lala
WORK_DIR=/tmpdir/wilson/dennis/$JOB_ID/
RESULTS_DIR=$HOME/results/$JOB_ID/

mkdir -p $RESULTS_DIR

i=0
for proc in $WORK_DIR/*
do
    cat $proc/log | grep 'S:' | rev | cut -d ':' -f 1 | rev > $RESULTS_DIR/s_$i.log
    cat $proc/log | grep 'R:' | rev | cut -d ':' -f 1 | rev > $RESULTS_DIR/rew_$i.log
    echo $i
    let i=i+1
done
