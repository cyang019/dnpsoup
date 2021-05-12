#!/bin/bash
#SBATCH -p sched_mit_psfc
#SBATCH --ntasks=102
#SBATCH --cpus-per-task=1
#SBATCH -J DNP
#SBATCH --export=OMP_NUM_THREADS,ALL
module rm gcc/4.8.4
module add gcc/8.3.0
module load gcc/8.3.0

# launch the code
export OMP_NUM_THREAD=1
export PATH="/home/cyang019/coding/projects/dnpsoup/build/dnpsoup_cli/:$PATH"
srun --multi-prog qband_p74_flip1_mix20k_d200k_loop200_scan_gB1_inc0p1.config
