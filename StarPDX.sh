#!/bin/bash

echo "#########"
echo "Adapter Trim"
echo 'Threads = ' $1
THR=$1

module purge
module load bbmap/37.28
ADAPTERS=/home/bgudenas/Annots/adapters.fa


R1=$(ls $TMPDIR/*R1.fastq )
ls $R1
R2=$(ls $TMPDIR/*R2.fastq )
ls $R2

SAMP=$( basename $R1 | rev | cut -c 10- | rev )

  bbduk.sh -Xmx32g threads=$THR in1=$R1 in2=$R2 \
        out1=$TMPDIR/${SAMP}_T1.fastq out2=$TMPDIR/${SAMP}_T2.fastq \
        ref=$ADAPTERS ktrim=r ordered \
        k=23 mink=11 hdist=1 tbo tpe \
        rcomp=f \
        ow=t \
        minlen=30 \
        trimq=5 \
        qtrim=rl
rm $TMPDIR/*_R*.fastq
[bgudenas@hpc02 PDX]$ cat ~/src/StarPDX.sh 
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


R1=$(ls $TMPDIR/*1.fastq )
ls $R1
R2=$(ls $TMPDIR/*2.fastq )
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

module load star/2.6.1c
GATK=/home/bgudenas/Software/gatk-4.1.1.0/gatk
module load java/1.8.0_181

if [ ! -f ./Data/Final/${SAMP}Aligned.out.bam ]
	then
$GATK CleanSam \
	--TMP_DIR=$TMPDIR \
	--VERBOSITY=ERROR \
	--VALIDATION_STRINGENCY=LENIENT \
	-I=$Human	\
	-O=$TMPDIR/${SAMP}.bam

$GATK SamToFastq \
    --TMP_DIR $TMPDIR \
    -F $TMPDIR/${SAMP}_R1.fastq \
    -F2 $TMPDIR/${SAMP}_R2.fastq \
    -I=$TMPDIR/${SAMP}.bam

ls -lh $R1
ls -lh $R2

  STAR --runMode alignReads --genomeDir $HsGENDIR --readFilesIn $TMPDIR/${SAMP}_R1.fastq $TMPDIR/${SAMP}_R2.fastq \
    --runThreadN $THR \
    --outSAMtype BAM  Unsorted \
    --quantMode GeneCounts \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.05 \
    --outFileNamePrefix ./Data/Final/${SAMP} \
    --outTmpDir $TMPDIR/tmp${SAMP} \
    --outFilterType BySJout \
    --peOverlapNbasesMin 10 \
    --outFilterMultimapNmax 20 \
    --alignSJDBoverhangMin 1

fi
