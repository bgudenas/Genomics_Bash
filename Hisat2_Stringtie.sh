#!/bin/bash

module purge
module load gcc
~/bin/hisat2 --version

echo "HISAT2 -- Stringtie"
mkdir -p Stringtie

echo 'Threads = ' $1
THR=$1
echo 'GENDIR = ' $2
GENDIR=$2
echo 'GTF = ' $3
GTF=$3
echo 'SAMP =' $4
SAMP=$4

ls -l $TMPDIR

OUTSAM=./Stringtie/${SAMP}.sam
SORTBAM=./Stringtie/${SAMP}_sort.bam
if [ ! -f $SORTBAM ];
then

if [ ! -f $OUTSAM ];
then
echo "HISAT2 Alignemnt"

R1=$(ls $TMPDIR/*1.fastq )
ls $R1
R2=$(ls $TMPDIR/*2.fastq )
ls $R2

hisat2 --dta -p $THR \
-x $GENDIR -1 $R1 -2 $R2 \
--rna-strandness RF \
-S $OUTSAM

fi

echo "SORTING"
module load samtools/1.9
samtools sort -@ $THR -m 40G -o $SORTBAM ./Stringtie/${SAMP}.sam
fi

OUT=./string_GTF/
mkdir -p $OUT

if [ ! -f $OUT/${SAMP}.gtf ]
then

stringtie -p $THR --rf -G $GTF -o $OUT/${SAMP}.gtf $SORTBAM
fi


if [ -f ./stringtie_merged.gtf ]
then
if [ ! -f ./Ballgown/${SAMP}.gtf ]
then

mkdir -p Ballgown/${SAMP}
stringtie -p $THR --rf -B -e -G ./stringtie_merged.gtf -o Ballgown/${SAMP}/${SAMP}.gtf $SORTBAM

	fi
fi
