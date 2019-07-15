#!/bin/bash

echo 'Threads = ' $1
THR=$1
echo 'SAMPNUM = ' $2
SAMPNUM=$2

module purge
module load star/2.6.1c
module load sambamba
module load samtools/1.7
GENDIR=/home/bgudenas/Annots/Human/star_GRCh38_99

R1=$(ls $TMPDIR/*1.fastq )
ls $R1
R2=$(ls $TMPDIR/*2.fastq )
ls $R2

SAMP=$( basename $R1 | rev | cut -c 10- | rev )

mkdir -p ./Stringtie/BAM
mkdir -p ./Stringtie/Samps
mkdir -p ./Stringtie/Merged
mkdir -p ./Stringtie/Ballgown

INBAM=$TMPDIR/${SAMP}Aligned.out.bam
SORTBAM=$TMPDIR/${SAMP}_Sorted.bam
XSBAM=./Stringtie/BAM/${SAMP}_XS.bam
if [ ! -f  $XSBAM ];
    then

STAR --runMode alignReads --genomeDir $GENDIR --readFilesIn $R1 $R2 \
    --runThreadN $THR \
    --outFileNamePrefix $TMPDIR/${SAMP} \
    --outSAMtype BAM Unsorted \
    --outSAMattrIHstart 0 \
    --twopassMode Basic \
    --alignSJstitchMismatchNmax 1 -1 1 1 \
    --outFilterType BySJout

sambamba sort -t $THR -m 38G --tmpdir $TMPDIR \
    -o $SORTBAM $INBAM

	samtools view --threads $THR -h $SORTBAM | awk -v strType=2 -f /home/bgudenas/src/tagXSstrandedData.awk | samtools view --threads $THR -bS - > $XSBAM
fi

	printf "STRINGTIE ASSEMBLY ############ \n"
	date
	stringtie --version
	GTF=/home/bgudenas/Annots/Human/Homo_sapiens.GRCh38.93.gtf
	OUT=./Stringtie/Samps/${SAMP}.gtf

if [ ! -f $OUT ];
	then
	stringtie -p $THR --rf -G $GTF -o $OUT -l $SAMP $XSBAM
fi


COUNT=$( ls ./Stringtie/Samps/*.gtf | wc -l )
if [ ! -f $MEGA ] && [ $SAMPNUM -eq $COUNT ] ;
	then
	MERGE="./Stringtie/Assembled/merge_list.txt"
	ls ./Stringtie/Samps/*.gtf > $MERGE
	cat $MERGE
	wc -l MERGE
	MEGA=./Stringtie/Merged/stringtie_merged.gtf

	stringtie --merge -p $THR -G $GTF -o $MEGA $MERGE 
fi

	GENE_ABS=./Stringtie/Ballgown/${SAMP}/gene_abundances.tsv
    OUT=./Stringtie/Ballgown/${SAMP}/${SAMP}.gtf
if [ -f $MEGA ] &&  [ ! -f $OUT];
	then
	stringtie -p $THR -B -e --rf -G $MEGA -A $GENE_ABS -o $OUT $XSBAM
fi
