#!/bin/bash

cd ~/rnafold

echo Job started running at `date` on `hostname`

export PYTHONPATH=${PYTHONPATH}:~/rnafold

echo Running at: `pwd`
echo Python path:
echo ${PYTHONPATH}
ls -l store_new_shuffles.py*

~/anaconda2/bin/dask-worker --scheduler-file ~/rnafold/dask-scheduler.json --nprocs 1 --nthreads 1 --reconnect --death-timeout 30

echo Job finished running at `date`

