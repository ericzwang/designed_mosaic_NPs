#!/bin/bash


for((i=0;i<=9;i++)); do
    cd $i
    sbatch --job-name="agg${i}" run.sh
    cd ..
done
