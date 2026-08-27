#!/bin/bash

#SBATCH --job-name mr
#SBATCH -q primary
#SBATCH --cpus-per-task=24
# Request memory
#SBATCH --mem-per-cpu=3GB
# Mail when the job begins, ends, fails, requeues
#SBATCH --mail-type=ALL
# Where to send email alerts
#SBATCH --mail-user=he0194@wayne.edu
# Create an output file that will be output_<jobid>.out
#SBATCH -o output_%j.out
# Create an error file that will be error_<jobid>.out
#SBATCH -e errors_%j.err
# Set maximum time limit
#SBATCH -t 30-0

# OpenMP re-enable (2026-08-27) is unproven at real scale and measured ~9%
# SLOWER than serial at small scale -- force serial for this run until the
# per-site ablation (see NOTES.md) settles which sites actually help.
export OMP_SITE_MASK=0

~/MinimalResolution/mr_st 270 240
~/MinimalResolution/BPtab 270
~/MinimalResolution/mr_BP 270 220

