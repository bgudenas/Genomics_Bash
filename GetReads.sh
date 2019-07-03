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
