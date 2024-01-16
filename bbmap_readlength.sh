#!/bin/bash
#PBS -l select=1:ncpus=6:mem=35gb:scratch_local=150gb
#PBS -l walltime=24:00:00
#PBS -N bbmap_readlength
#
#Creates histogram of read length

INPUT=/full/path/to/Physko22092023_bonito_v073/bonito_v073_KO22.fq
OUTFILE=./bonito_v073_KO22_histogram.txt
OUTDIR=/full/path/to/Physko22092023_bonito_v073/

cp $INPUT $SCRATCH/
INPUT=$SCRATCH/$(basename $INPUT)
cd $SCRATCH


source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module load bbmap


readlength.sh in=$INPUT out=$OUTFILE
stats.sh in=$INPUT

mkdir -p $OUTDIR
cp $OUTFILE $OUTDIR/
rm -r $SCRATCH/*