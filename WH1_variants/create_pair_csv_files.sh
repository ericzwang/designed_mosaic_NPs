#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1


source /home/gridsan/ezw/.bashrc
module load anaconda

python create_pair_csv_files.py $1 >& logfiles/log$1.out

