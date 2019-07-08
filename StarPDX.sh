#!/bin/bash

echo "starPDX"
THR=$1

module purge
module load star/2.6.1c
module load sambamba

HsGENDIR="/home/bgudenas/Annots/Human/star_GRCh38_99/"
MmGENDIR="/home/bgudenas/Annots/Mouse/mm10_star99/"
mkdir -p ./Data/Human
mkdir -p ./Data/Mouse


R1=$(ls $TMPDIR/*T1.fastq )
ls $R1
R2=$(ls $TMPDIR/*T2.fastq )
ls $R2

SAMP=$( basename $R1 | rev | cut -c 10- | rev )


### Human
INBAM="./Data/Human/${SAMP}Aligned.out.bam"
SORTBAM="./Data/Human/${SAMP}_Sorted.bam"

if [ ! -f  $SORTBAM ]; ## if Bam does not exist
    then
    echo "STAR"
    STAR --runMode alignReads --genomeDir $HsGENDIR --readFilesIn $R1 $R2 \
    --runThreadN $THR \
    --outSAMtype BAM  Unsorted \
    --outSAMattributes NM \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.05 \
    --outFileNamePrefix ./Data/Human/${SAMP} \
    --outTmpDir $TMPDIR/tmp${SAMP} \
    --outFilterType BySJout \
    --peOverlapNbasesMin 10 \
    --outFilterMultimapNmax 20 \
    --alignSJDBoverhangMin 1

    sambamba sort -t $THR -m 38G --tmpdir $TMPDIR \
	--sort-by-name \
	-o $SORTBAM $INBAM 
#  rm $INBAM
fi

##Mouse
INBAM="./Data/Mouse/${SAMP}Aligned.out.bam"
SORTBAM="./Data/Mouse/${SAMP}_Sorted.bam"
if [ ! -f  $SORTBAM ]; ## if Bam does not exist
    then
    echo "STAR"
    STAR --runMode alignReads --genomeDir $MmGENDIR --readFilesIn $R1 $R2 \
    --runThreadN $THR \
    --outSAMtype BAM  Unsorted \
    --outSAMattributes NM \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.05 \
    --outFileNamePrefix ./Data/Mouse/${SAMP} \
    --outTmpDir $TMPDIR/tmp${SAMP} \
    --outFilterType BySJout \
    --peOverlapNbasesMin 10 \
    --outFilterMultimapNmax 20 \
    --alignSJDBoverhangMin 1

    sambamba sort -t $THR -m 38G --tmpdir $TMPDIR \
    --sort-by-name \
    -o $SORTBAM $INBAM
 # rm $INBAM
fi

Human=./Disambig/${SAMP}.disambiguatedSpeciesA.bam
if [ ! -f $Human ]
	then
module load disambiguate
disambiguate.py -s ${SAMP} \
		./Data/Human/${SAMP}_Sorted.bam \
		./Data/Mouse/${SAMP}_Sorted.bam  \
		-i $TMPDIR \
		-a star \
		-d \
		-o Disambig
fi

module purge
module load htseq
which python

GTF=/home/bgudenas/Annots/Human/Homo_sapiens.GRCh38.93.gtf
htseq-count -f bam -r name -s yes \
$Human $GTF > Counts/${SAMP}_Counts.txt
