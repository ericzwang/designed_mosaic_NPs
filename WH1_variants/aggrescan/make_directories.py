import os
import pandas as pd

distance = 5
num_jobs = 300
size_per_job = int(pd.read_csv('../csv_files/stablegroup1.csv').shape[0]/num_jobs) + 1
slurm_file_string = """#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1

source /home/gridsan/ezw/.bashrc
module load anaconda
source activate aggrescan


"""

def make_WT_dir():
    os.mkdir('jobs/WT')
    with open('jobs/WT/scores.csv', 'w') as f:
        f.write('')
    with open('jobs/WT/run.sh', 'w') as f:
        f.write(slurm_file_string)
        f.write(f'aggrescan -i ../../RBD.pdb --chain C -f ~/foldx -D {distance} -v 3 -r \n\n')
        f.write("awk -F ',' '{print $5}' A3D.csv | tail -n +2 | awk '{ sum += $1 } END { print sum }' | bc >> scores.csv\n\n")


def format_muts(jobn, group, num_jobs=300):
    df_group = pd.read_csv(f'../csv_files/stablegroup{group}.csv')
    muts_strs = []
    inds = list(range(jobn*size_per_job, (jobn+1)*size_per_job))
    if jobn == num_jobs - 1:
        inds = list(range(jobn*size_per_job, df_group.shape[0]))
    for ind in inds:
        muts_reformat = [mut[0]+mut[-1]+mut[1:-1]+'C' for mut in df_group.iloc[ind]]
        muts_str = " ".join(muts_reformat)
        muts_strs.append(muts_str)
    return muts_strs

def make_job_dir(jobn):
    os.mkdir('jobs/'+str(jobn))
    with open(f'jobs/{jobn}/scores.csv', 'w') as f:
        f.write('')
    with open(f'jobs/{jobn}/run.sh', 'w') as f:
        f.write(slurm_file_string)
        for group in [1,2]:
            muts = format_muts(jobn,group)
            for i,seq in enumerate(muts):
                ind = jobn*size_per_job + i
                f.write(f'aggrescan --mutate {seq} -i ../../RBD.pdb --chain C -f ~/foldx -D {distance} -v 3 -r \n')
                f.write("score=`awk -F ',' '{print $5}' A3D.csv | tail -n +2 | awk '{ sum += $1 } END { print sum }' | bc`\n")
                f.write(f"echo {ind}_{group},$score >> scores.csv\n\n")

def run():
    make_WT_dir()
    for jobn in range(num_jobs):
        make_job_dir(jobn)



if __name__ == '__main__':
    run()