#!/bin/bash

echo "#########"
echo "Kallisto"
echo 'Threads = ' $1
THR=$1

module purge
module load kallisto/0.43.1 
GENDIR="/home/bgudenas/Annots/Human/Kallisto/GRCh38_transcripts.idx"
mkdir -p Kallisto


R1=$(ls $TMPDIR/*1.fastq )
ls $R1
R2=$(ls $TMPDIR/*2.fastq )
ls $R2

SAMP=$( basename $R1 | rev | cut -c 10- | rev )

if [ ! -f ./Kallisto/${SAMP}/abundance.tsv ];
	then
kallisto quant -t $THR -i $GENDIR --bias --rf-stranded -o ./Kallisto/${SAMP} -b 100 $R1 $R2

fi
