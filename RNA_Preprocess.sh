#!/bin/bash

echo "#########"
echo "Adapter Trim"
echo 'Threads = ' $1
THR=$1

module purge
module load bbmap/37.28
ADAPTERS=/home/bgudenas/Annots/adapters.fa

R1=$(ls $TMPDIR/*T1.fastq )
ls $R1
R2=$(ls $TMPDIR/*T2.fastq )
ls $R2

SAMP=$( basename $R1 | rev | cut -c 10- | rev )


GENDIR=/research/rgs01/project_space/northcgrp/Northcott_Bioinformatics/northcgrp/References/Human/GATK/hg38/

STAR --runMode alignReads --genomeDir $GENDIR --sjdbOverhang 99 --readFilesIn $R1 $R2 \
    --runThreadN $THR \
    --outFileNamePrefix ./Data/${SAMP} \
    --outSAMtype BAM Unsorted \
    --outFilterType BySJout \
    --twopassMode Basic \
    --limitSjdbInsertNsj 2000000 \
    --limitOutSJcollapsed 2000000 \
    --outTmpDir $TMPDIR/tmp${SAMP}


     sambamba sort -t $THR -m 38G --tmpdir $TMPDIR \
    --sort-by-name \
    -o $SORTBAM $INBAM



java -jar -XX:ParallelGCThreads=1 \
	$PICARD FixMateInformation \
	I=$BAM \
	O=./GATK/${BASE}/${SAMP}_MC.bam \
	ADD_MATE_CIGAR=true \
	SO=coordinate \
	VALIDATION_STRINGENCY=LENIENT \
	TMP_DIR=$TMPDIR

java -jar -Djava.io.tmpdir=$TMPDIR -XX:ParallelGCThreads=1 \
		$PICARD AddOrReplaceReadGroups \
		I=./GATK/${BASE}/${SAMP}_region.bam \
		O=./GATK/${BASE}/${SAMP}_RG.bam \
		VALIDATION_STRINGENCY=LENIENT \
		SO=coordinate \
		RGID=$SAMP RGSM=$SAMP RGLB=$SAMP RGPL=$PLATFORM RGPU=$MACHINE

		fi

	if [ ! -f ./GATK/${BASE}/${SAMP}_dup.bam ]
	 	then

java -jar -Djava.io.tmpdir=$TMPDIR -XX:ParallelGCThreads=1 \
		$PICARD MarkDuplicates \
		I=./GATK/${BASE}/${SAMP}_RG.bam \
		O=./GATK/${BASE}/${SAMP}_dup.bam \
		M=./GATK/${BASE}/${SAMP}_metric.txt \
		CREATE_INDEX=true \
		VALIDATION_STRINGENCY=LENIENT

		fi

if [ ! -f ./GATK/${BASE}/${SAMP}_split.bam ]
	 	then

java -jar -Djava.io.tmpdir=$TMPDIR -XX:ParallelGCThreads=1 \
		$GATK -T SplitNCigarReads \
		-R $FASTA \
		-I ./GATK/${BASE}/${SAMP}_dup.bam \
		-o ./GATK/${BASE}/${SAMP}_split.bam \
		-rf ReassignOneMappingQuality \
		-RMQF 255 -RMQT 60  \
		-U ALLOW_N_CIGAR_READS


java -jar -Djava.io.tmpdir=$TMPDIR -XX:ParallelGCThreads=1 \
		$GATK -T BaseRecalibrator \
			-R $FASTA \
			-I ./GATK/${BASE}/${SAMP}_split.bam \
			-o ./GATK/${BASE}/${SAMP}_recal.table

java -jar -Djava.io.tmpdir=$TMPDIR -XX:ParallelGCThreads=1 \
		$GATK -T PrintReads \
			-R $FASTA \
			-BQSR ./GATK/${BASE}/${SAMP}_recal.table \
			-I ./GATK/${BASE}/${SAMP}_split.bam \
			-o ./GATK/${BASE}/${SAMP}_BQSR.bam
