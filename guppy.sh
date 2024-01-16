#!/bin/bash -ex
#PBS -q gpu@meta-pbs.metacentrum.cz
#PBS -l walltime=7:00:00
#PBS -l select=1:ncpus=1:ngpus=1:mem=30gb:gpu_mem=20gb:scratch_local=100gb 
#PBS -N guppy_6_3_8


DATADIR=/full/path/to/gupy_report
mkdir -p $DATADIR

cd $SCRATCH

echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" >> $DATADIR/jobs_info.txt

source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module add guppy-6.3.8-gpu

cp -r /full/path/to/fast5_pass/ ./ || { echo >&2 "Error while copying input file(s)!"; exit 2; }

guppy_basecaller -i ./fast5_pass -r -s ./guppy_v6_3_8_KO22_davka -c dna_r10.4.1_e8.2_400bps_hac.cfg -x auto --gpu_runners_per_device 16 --num_callers 16 --chunks_per_runner 2000 --trim_strategy none --disable_qscore_filtering

mv ./guppy_v6_3_8_KO22_davka $DATADIR/

rm -r $SCRATCH/*