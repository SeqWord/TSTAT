#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=30:00:00
#PBS -q normal
#PBS -o path
#PBS -e path
#PBS -k oe
#PBS -m ae
#PBS -M your@e-mail
#PBS -N tstat

module load python-3.9.5
module load R-4.4.3

python3 /path/TSTAT/server_run.py \
    --project test \        # Project directory within 'input' directory
    --embedded '' \         # Sub-project directory or ''
    --reference test.gbk \  # Reference GBK file
    --type single \         # DNA read type: single | paired-end | smart
    --basename K8B78 \      # Locus tag specification in the reference GBK filr
    --control C \           # Unique marker of control fastq files in their names
    --experiment E \        # Unique marker of experiment fastq files in their names
    --plots "Volcano plot,BaseMean plot"
