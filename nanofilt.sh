#!/bin/bash
#PBS -l select=1:ncpus=6:mem=35gb:scratch_local=150gb
#PBS -l walltime=10:00:00
#PBS -N nanofilt_2
#

INPUT=/full/path/to/Physko22092023_bonito_v703/bonito_v073_KO22.fq
OUTDIR=/full/path/tp/katka/Physko22092023_bonito_v703/
OUTFILE=./all_bonito_v073_KO22_trimmed.fastq

cd $SCRATCH
cp $INPUT $SCRATCH/
INPUT=$(basename $INPUT)



source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module add python36-modules-gcc

NanoFilt -l 500 --headcrop 10 < $INPUT > $OUTFILE

mkdir -p $OUTDIR
cp $OUTFILE $OUTDIR/



rm -r $SCRATCH/*