#!/bin/bash

echo 'Threads = ' $1
THR=$1
echo 'SAMPNUM = ' $2
SAMPNUM=$2
echo 'GENDIR = ' $3
GENDIR=$3

module purge
module load star/star/2.7.1a
module load sambamba

GENDIR=/home/bgudenas/Annots/Human/star_GRCh38_99

R1=$(ls $TMPDIR/*1.fastq )
ls $R1
R2=$(ls $TMPDIR/*2.fastq )
ls $R2

SAMP=$( basename $R1 | rev | cut -c 10- | rev )
echo $SAMP

mkdir -p ./Quant

#INBAM=$TMPDIR/${SAMP}Aligned.out.bam
#SORTBAM=$TMPDIR/${SAMP}_SO.bam
if [ ! -f  ./Quant/${SAMP}Log.final.out ];
    then

STAR --runMode alignReads --genomeDir $GENDIR --readFilesIn $R1 $R2 \
    --runThreadN $THR \
    --outFileNamePrefix $TMPDIR/${SAMP} \
    --outSAMtype BAM Unsorted \
    --outFilterType BySJout \
    --outFilterMismatchNmax 999 \
	--outFilterMismatchNoverReadLmax 0.05 \
	--outFilterMultimapNmax 20 \
	--alignSJDBoverhangMin 1 \
	--outTmpDir $TMPDIR/tmp${SAMP}

mv $TMPDIR/${SAMP}*Reads* ./Quant
mv $TMPDIR/${SAMP}*Log.final* ./Quant

#sambamba sort -t $THR -m 38G --tmpdir $TMPDIR \
#-o $SORTBAM $INBAM

module purge
module load kallisto/0.43.1 
GENDIR="/home/bgudenas/Annots/Human/Kallisto/GRCh38_transcripts.idx"
mkdir -p Kallisto


if [ ! -f ./Kallisto/${SAMP}/abundance.tsv ];
	then
kallisto quant -t $THR -i $GENDIR --bias --rf-stranded -o ./Kallisto/${SAMP} -b 100 $R1 $R2

fi

