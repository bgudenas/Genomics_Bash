#!/bin/bash

echo "#########"
echo "Get Reads"
echo 'FDIR = ' $1
FDIR=$1
echo 'NCUT = ' $2
NCUT=$2
echo 'HCUT = ' $3
HCUT=$3
echo 'LSB_JOBINDEX = ' $4
LSB_JOBINDEX=$4
echo 'RUN = ' $5
RUN=$5

FASTQS=$( find ${FDIR}* -type f -name "*q.gz" -exec ls {} + )
SAMP=$( ls $FASTQS | rev | cut -c ${NCUT}- | rev  |sed 's!.*/!!' | sed 's/_$//' | sort | uniq -d | awk -v var="$LSB_JOBINDEX" 'FNR == var {print}' )

[[ !  -z  ${SAMP// /}  ]] || { echo >&2 "SAMPLE ID is EMPTY"; exit 1; }

## Cut sample run ID header if needed
if (( HCUT > 0 ));
then
SAMP=$( echo $SAMP | cut -c ${HCUT}- )
fi

printf "\n"
echo "SAMPLE####################"
echo $SAMP
R1=$( ls $FASTQS| grep ${SAMP} | grep "R1" )
R2=$(  ls $FASTQS  | grep ${SAMP} | grep "R2" )

## Abort if not same number of read1 and read2 fastqs
if (( $( ls $R1 | wc -l ) != $( ls $R2 | wc -l ) )); 
then
echo >&2 "FASTQS DONT MATCH"
exit 1
fi

echo "READ1####"
ls -lh $R1
echo "READ2####"
ls -lh $R2

MR1="${TMPDIR}/${SAMP}_R1.fastq"
MR2="${TMPDIR}/${SAMP}_R2.fastq"

if [ "$RUN" = "YES" ]
	then
cat $R1 | gunzip > $MR1
cat $R2 | gunzip > $MR2
fi
[bgudenas@hpc01 PDX]$ cat ~/src/StarPDX.sh 
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

    sambamba sort -t $THR -m 45G --tmpdir $TMPDIR \
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

    sambamba sort -t $THR -m 45G --tmpdir $TMPDIR \
	--sort-by-name \  
	-o $SORTBAM $INBAM 
 # rm $INBAM
fi


module load disambiguate
disambiguate.py -s ${SAMP} \
		./Data/Human/${SAMP}_Sorted.bam \
		./Data/Human/${SAMP}_Sorted.bam  \
		-i $TMPDIR \
		-a star \
		-d \
		-o Disambig
