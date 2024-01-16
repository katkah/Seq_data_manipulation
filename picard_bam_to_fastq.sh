#!/bin/bash
#PBS -l select=1:ncpus=6:mem=35gb:scratch_local=150gb
#PBS -l walltime=10:00:00
#PBS -N picard_telomeres
#
#Converts bam to fastq


INPUT=/full/path/to/Physko22092023_bonito_v073/mapping_reads_with_telomere_motives/ko22_bonito_v073_map_F4_sort.bam
OUTDIR=/full/path/to/Physko22092023_bonito_v073/mapping_reads_with_telomere_motives/
OUTFILE=./ko22_bonito_v073_map_F4_sort.fq


cp $INPUT $SCRATCH/
INPUT=$SCRATCH/$(basename $INPUT)
cd $SCRATCH/


cd $SCRATCH

source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module load picard
module load jdk


java -jar /software/picard/2.9.0/build/libs/picard.jar SamToFastq I=$INPUT F=$OUTFILE


mkdir -p $OUTDIR
cp $OUTFILE $OUTDIR/
rm -r $SCRATCH/*