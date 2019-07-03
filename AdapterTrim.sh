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
