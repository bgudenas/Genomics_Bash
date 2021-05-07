#!/bin/bash

echo "#########"
echo "Get Reads"
echo 'FDIR = ' $1
FDIR=$1

echo 'SAMP = ' $2
SAMP=$2

FASTQS=$( find ${FDIR}* -type f -name "*q.gz" -exec ls {} + )

##SAMP=$(  ls $FASTQS |  rev | cut -c ${NCUT}- | rev  | sed 's!.*/!!' | sed 's/_$//' | sort | uniq -d | awk -v var="$LSB_JOBINDEX" 'FNR == var {print}' )

####"
echo $SAMP
R1=$( ls $FASTQS| grep ${SAMP} | grep "_R1" )
R2=$(  ls $FASTQS  | grep ${SAMP} | grep "_R2" )

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

cat $R1 | gunzip > $MR1
cat $R2 | gunzip > $MR2
