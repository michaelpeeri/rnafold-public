#!/bin/bash

cd /tamir1/mich1/
source .bashrc

cd ~/cellfold

echo Job started running at `date` on `hostname`

export PYTHONPATH=${PYTHONPATH}:~/cellfold

echo Running at: `pwd`
echo Python path:
echo ${PYTHONPATH}

~/anaconda3/bin/dask-worker --scheduler-file ~/cellfold/dask-scheduler.json --nprocs 1 --nthreads 1 --reconnect --death-timeout 30

echo Job finished running at `date`

