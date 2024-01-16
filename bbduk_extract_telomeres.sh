#PBS -l walltime=24:00:00
#PBS -N bbmap
#
#Creates a file containing telomere reads

INPUT=/full/path/to/Physko22092023_bonito_v703/bonito_v073_KO22_trim_map_F4_sort.fq
OUTFILE=./bonito_v073_KO22_trim_map_F4_sort_tel.fq
OUTDIR=/full/path/to/Physko22092023_bonito_v703/

cp $INPUT $SCRATCH/
INPUT=$SCRATCH/$(basename $INPUT)
cd $SCRATCH


source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module load bbmap

#bbduk.sh in=YOUR.FASTQ outm=MATCHED.FASTQ \
#literal=SEARCH_STRING k=STRING_LENGTH hdist=NUMBER_OF_ERRORS

#You can specify whether or not BBDuk looks for the reverse-complement of 
#the reference sequences as well as the forward sequence with the flag “rcomp=t” or “rcomp=f”;
# by default it looks for both.

bbduk.sh in=$INPUT outm=$OUTFILE literal=TTTAGGGTTTAGGGTTTAGGG k=21 hdist=1

mkdir -p $OUTPUTDIR
cp $OUTFILE $OUTPUTDIR/