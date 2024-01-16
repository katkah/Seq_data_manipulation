#!/bin/bash -ex
#PBS -q gpu@meta-pbs.metacentrum.cz
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=1:ngpus=1:mem=30gb:gpu_mem=20gb:scratch_local=100gb 
#PBS -N bonito_v0_7_3

DATADIR=/full/path/to/Physko22092023_bonito_v703
mkdir -p $DATADIR
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" >> $DATADIR/jobs_info.txt
INPUTDIR=/full/path/to/fast5_pass/


cd $SCRATCH
mkdir $SCRATCH/reads
cp $INPUTDIR/* ./reads || { echo >&2 "Error while copying input file(s)!"; exit 2; }


module add mambaforge
mamba activate /full/path/to/home/user/bonito/
bonito basecaller dna_r10.4.1_e8.2_400bps_hac@v4.2.0 $SCRATCH/reads > bonito_v703_KO22.bam

mv ./bonito_v703_KO22.bam $DATADIR/

rm -r $SCRATCH/*