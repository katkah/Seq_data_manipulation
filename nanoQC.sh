#!/bin/bash
#PBS -l select=1:ncpus=6:mem=35gb:scratch_local=150gb
#PBS -l walltime=10:00:00
#PBS -N nanoQC
#

INPUT=/full/path/to/Physko22092023_bonito_v703/bonito_v073_KO22.fq
OUTPUT_DIR=/full/path/to/Physko22092023_bonito_v703/nanoQC_basecalled_data


cp $INPUT $SCRATCH/
INPUT=$SCRATCH/$(basename $INPUT)
cd $SCRATCH/

mkdir -p $SCRATCH/nanoQC

module add python36-modules-gcc

nanoQC -o $SCRATCH/nanoQC $INPUT

cd $SCRATCH/nanoQC

mkdir -p $OUTPUT_DIR

cp -r $SCRATCH/nanoQC $OUTPUT_DIR/

rm -r $SCRATCH/*