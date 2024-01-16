#!/bin/bash
#PBS -l select=1:ncpus=6:mem=35gb:scratch_local=150gb
#PBS -l walltime=10:00:00
#PBS -N minimap2_telomeres
#
# Runs minimap2

INPUT=/full/path/to/Physko22092023_guppy/guppy_v638_KO22_telomeres.fq
GENOME=/full/path/to/nanopore_basecalling/Physco_genome/Phytozome/PhytozomeV10/early_release/Ppatens_318_v3.3/assembly/Ppatens_318_v3.fa.gz
OUTDIR=/full/path/to/Physko22092023_guppy/mapping_reads_with_telomere_motives
OUTFILE=./guppy_v638_KO22_telomeres.sam
cp $INPUT $SCRATCH/
cp $GENOME $SCRATCH/


cd $SCRATCH


source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module load minimap2
echo "`minimap2 --version`"

minimap2 -ax map-ont $GENOME $INPUT > $OUTFILE

mkdir -p $OUTDIR
cp $OUTFILE $OUTDIR/

module unload minimap2
source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module load samtools

samtools flagstat $OUTFILE

#############################
#extract only mapped reads

INPUT=$OUTFILE
OUTFILE=./guppy_v638_KO22_telomeres_F4.sam



samtools view -h -F 4 $INPUT > $OUTFILE 

cp $OUTFILE $OUTDIR/

################################
#SORT AND INDEX
INPUT=$OUTFILE
OUTFILE=./guppy_v638_KO22_telomeres_F4_sort.bam

samtools sort $INPUT -o $OUTFILE
samtools index $OUTFILE


cp $OUTFILE $OUTDIR/
cp *.bai $OUTDIR/

#CLEAN SCRATCH
rm -r $SCRATCH/*