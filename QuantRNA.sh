#!/bin/bash
# TODO add option to retain sorted bam file

echo 'Threads = ' $1
THR=$1
echo 'GENDIR = ' $2
GENDIR=$2

module purge
module load star/2.7.1a
module load sambamba

echo "TMPDIR"
ls $TMPDIR
R1=$(ls $TMPDIR/*1.fastq )
ls $R1
R2=$(ls $TMPDIR/*2.fastq )
ls $R2

SAMP=$( basename $R1 | rev | cut -c 10- | rev )
echo $SAMP

mkdir -p ./Quant

INBAM=$TMPDIR/${SAMP}Aligned.out.bam
SORTBAM=./Quant/${SAMP}_SO.bam

if [ ! -f  ./Quant/${SAMP}Log.final.out ];
    then

STAR --runMode alignReads --genomeDir $GENDIR --readFilesIn $R1 $R2 \
	--runThreadN $THR \
	--outFileNamePrefix $TMPDIR/${SAMP} \
	--outSAMtype BAM Unsorted \
	--quantMode GeneCounts \
	--outFilterType BySJout \
	--outFilterMismatchNmax 999 \
	--outFilterMismatchNoverReadLmax 0.05 \
	--outFilterMultimapNmax 20 \
	--alignSJDBoverhangMin 1 \
	--outTmpDir $TMPDIR/tmp${SAMP}

ls $TMPDIR
mv $TMPDIR/${SAMP}*Reads* ./Quant
mv $TMPDIR/${SAMP}*Log.final* ./Quant

#module load sambamba/0.7.1 
#sambamba sort -t $THR -m 38G --tmpdir $TMPDIR -o $SORTBAM $INBAM

fi

module purge
module load kallisto/0.43.1 
mkdir -p Kallisto
if [[ "$GENDIR" == *"Mouse"* ]];
then
	echo "Switching to Kallisto Mouse Index"
	GENDIR="/home/bgudenas/Annots/Mouse/Kallisto/Mm10_transcripts.idx"
elif [[ "$GENDIR" == *"Human"* ]];
then
	echo "Switching to Kallisto Human Index"
	GENDIR="/home/bgudenas/Annots/Human/Kallisto/GRCh38_transcripts.idx"
else
	echo "No index found"
fi

if [ ! -f ./Kallisto/${SAMP}/abundance.tsv ];
	then
kallisto quant -t $THR -i $GENDIR --bias --rf-stranded -o ./Kallisto/${SAMP} -b 100 $R1 $R2

fi
