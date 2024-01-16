#!/bin/bash
#PBS -l select=1:ncpus=6:mem=35gb:scratch_local=150gb
#PBS -l walltime=24:00:00
#PBS -N makeblastdb
#

INPUT=/full/path/to/Physko22092023_guppy/all_guppy_v6_3_8.fastq
QUERY=/full/path/to/hygromycin_blast/query.fasta
OUTDIR=/full/path/to/hygromycin_blast/

cp $INPUT $SCRATCH/
INPUT=$SCRATCH/$(basename $INPUT)
cp $QUERY $SCRATCH/
QUERY=$SCRATCH/$(basename $QUERY)

cd $SCRATCH

source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module load seqtk
#fastq to fasta
seqtk seq -a $INPUT > ${INPUT%.fastq}.fasta
module unload seqtk



module load blast-plus/2.12.0-gcc-8.3.0-ohlv7t4

makeblastdb -in ${INPUT%.fastq}.fasta -out all_guppy_db -dbtype 'nucl' -hash_index -max_file_sz '20GB'

blastn -query $QUERY -task blastn -db all_guppy_db -out all_guppy_blast.html -evalue 10 -word_size 4 -num_threads 6 -html

rm $INPUT
cp -r * $OUTDIR/