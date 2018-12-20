#BSUB -P SKBTB
#BSUB -R "rusage[mem=5000]"
#BSUB -n 12
#BSUB -oo ./Logs/KBTBD4_STAR.log -eo ./Logs/KBTBD4_STAR.err
#BSUB -R "span[hosts=1]"

RUN="YES"
NCUT=26
BASE="north"
TRIM="NO"
REPAIR="NO"
SPECIES="MOUSE"
PROJ="KBTBD4"
THR=12

## Map reads using STAR multi-sample protocol
## first map using basic params
## Store SJ junctions
## Second pass using all sample junctions

### fastq files should be in ./Raw/${BASE}/*/ with sample *.fastq.gz inside

date
module purge
module load bbmap/37.28
module load samtools/1.7
module load sambamba/0.5.5
printf "STAR RNA-seq mapping############# \n"
STAR --version


ADAPTERS=/home/bgudenas/Annots/truseq.fa

### ENSEMBL annots
GENDIR="/home/bgudenas/Annots/Mouse/mm10_star99"
if [ "$SPECIES" = "HUMAN" ]
    then
GENDIR="/home/bgudenas/Annots/Human/star_GRCh38_99"
fi
#########################################

SAMPNUM=$(ls ./Raw/${BASE}*/*/*.fastq.gz | rev | cut -c ${NCUT}- | rev  |  sed 's!.*/!!' | sort | uniq -d | wc -l)
printf "SAMPLES DETECTED = ${SAMPNUM} \n"
printf "FASTQ'S DETECTED = $(ls -l ./Raw/${BASE}*/*/* | wc -l) \n"

RLENGTH=$(gunzip -c "$(ls ./Raw/${BASE}*/*/*.gz | head -n1)" | head -n 1000 | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c)
printf "READ LENGTH(n=250) = ${RLENGTH} \n\n"

mkdir -p ./Data/${BASE}
mkdir -p ./Data/${BASE}/SJ
mkdir -p ./BAM/${BASE}


for i in $(ls ./Raw/${BASE}*/*/*.fastq.gz | rev | cut -c ${NCUT}- | rev  |sed 's!.*/!!' | sort | uniq )
do	
    printf "\n"  
    SAMP=$(basename ${i})
    echo "SAMPLE##################"
    echo $SAMP
    R1="./Raw/${BASE}*/*/"${SAMP}"*R1*.fastq.gz"
    R2="./Raw/${BASE}*/*/"${SAMP}"*R2*.fastq.gz"
    echo "READ1####"
    ls -lh $R1
    echo "READ2####"
    ls -lh $R2

###################################
    MR1="${TMPDIR}/Read1_${SAMP}.fastq.gz"
    MR2="${TMPDIR}/Read2_${SAMP}.fastq.gz"
    
if [ "$RUN" = "YES" ]
	then
    cat $R1 > $MR1
    cat $R2 > $MR2
    
    ls -lh $MR1
    ls -lh $MR2

if [ "$TRIM" = "YES" ]
    then

TR1="${TMPDIR}/Read1_Trim_${SAMP}.fastq.gz"
TR2="${TMPDIR}/Read2_Trim_${SAMP}.fastq.gz"

if [ "$REPAIR" = "YES" ]
    then
## Fix misordered paired-end files
FR1="${TMPDIR}/Read1_Fix_${SAMP}.fastq.gz"
FR2="${TMPDIR}/Read2_Fix_${SAMP}.fastq.gz"

repair.sh -Xmx12g in1=$MR1 in2=$MR2 out1=$FR1 out2=$FR2 outs=$TMPDIR/singletons.fq

## overwrite input to trim
MR1="${TMPDIR}/Read1_Fix_${SAMP}.fastq.gz"
MR2="${TMPDIR}/Read2_Fix_${SAMP}.fastq.gz"
ls -lh $MR1
ls -lh $MR2
fi

  bbduk.sh -Xmx12g in1=$MR1 in2=$MR2 \
      	out1=$TR1 out2=$TR2 \
      	ref=$ADAPTERS ktrim=r ordered \
       	k=23 mink=11 hdist=1 tbo tpe \
       	rcomp=f \
	    ow=t \
	    minlen=30 \
	    trimq=10 \
	    qtrim=rl

echo "TRIMMED ##########"
MR1="${TMPDIR}/Read1_Trim_${SAMP}.fastq.gz"
MR2="${TMPDIR}/Read2_Trim_${SAMP}.fastq.gz"
ls -lh $MR1
ls -lh $MR2
fi

if [ ! -f ./Data/${BASE}/SJ/${SAMP}SJ.out.tab ];
	then
    echo "STAR"
    STAR --runMode alignReads --genomeDir $GENDIR --readFilesIn $MR1 $MR2 \
    --runThreadN $THR --readFilesCommand zcat \
    --outFileNamePrefix $TMPDIR/${SAMP} \
    --outTmpDir $TMPDIR/tmp${SAMP} \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.05 \
    --alignIntronMin 20 \
    --alignSJDBoverhangMin 1 \
    --outFilterType BySJout \
    --peOverlapNbasesMin 10

    mv $TMPDIR/*SJ.out* ./Data/${BASE}/SJ
	fi
fi
done

##########################################################################################################
## Second pass using all indexes --sjdbFileChrStartEnd
## Move junction files to prevent overwriting
if [ "$RUN" = "YES" ]
	then

awk -f /home/bgudenas/src/sjCollapseSamples.awk ./Data/${BASE}/SJ/*SJ.out.tab | sort -k1,1V -k2,2n -k3,3n > ./Data/${BASE}/SJ/SJ.all
grep -E -v "KI|GL|MT" ./Data/${BASE}/SJ/SJ.all > ./Data/${BASE}/SJ/SJ.pass1
awk '{ if ($10 >= 2) { print } }' ./Data/${BASE}/SJ/SJ.pass1 > ./Data/${BASE}/SJ/SJ.Filtered

	fi

for i in $(ls ./Raw/${BASE}*/*/*.fastq.gz | rev | cut -c ${NCUT}- | rev  |sed 's!.*/!!' | sort | uniq -d)
do
    SAMP=$(basename ${i})
    echo "ROUND 2 -- ${SAMP}"
if [ "$RUN" = "YES" ]
	then
    MR1="${TMPDIR}/Read1_${SAMP}.fastq.gz"
    MR2="${TMPDIR}/Read2_${SAMP}.fastq.gz"
    if [ "$TRIM" = "YES" ]
        then
        MR1="${TMPDIR}/Read1_Trim_${SAMP}.fastq.gz"
        MR2="${TMPDIR}/Read2_Trim_${SAMP}.fastq.gz"
        fi

STAR --runMode alignReads --genomeDir $GENDIR --readFilesIn $MR1 $MR2 \
    --runThreadN $THR --readFilesCommand zcat \
    --outFileNamePrefix ./Data/${BASE}/${SAMP} \
    --outSAMtype BAM Unsorted \
    --quantMode GeneCounts  \
    --sjdbFileChrStartEnd ./Data/${BASE}/SJ/SJ.Filtered ./Data/${BASE}/SJ/${SAMP}SJ.out.tab \
    --peOverlapNbasesMin 10 \
    --alignSJDBoverhangMin 1 \
    --limitSjdbInsertNsj 2000000 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.05 \
    --alignIntronMin 20 \
    --outFilterType BySJout \
    --outTmpDir $TMPDIR/tmp2${SAMP}

   INBAM="./Data/${BASE}/${SAMP}Aligned.out.bam"
   SORTBAM="./BAM/${BASE}/${SAMP}Aligned.sortedByCoord.out.bam"

sambamba sort -t $THR -m 48G --tmpdir $TMPDIR \
 	-o $SORTBAM $INBAM   
	fi
done

cat ./Data/${BASE}/*Log.final* | grep -E 'Number of input reads|Uniquely mapped reads %'
date


if [ "$RUN" = "YES" ]
	then
	module load R/3.4.0
	cd ./Data/${BASE}
	Rscript /home/bgudenas/src/Parse_STARLogs.R ${BASE}
	Rscript /home/bgudenas/src/countMat.R "Human"
    cd ../..
    R -e "rmarkdown::render('/home/bgudenas/src/counts_EDA.Rmd', output_file = '${BASE}_EDA.html', output_dir = './Results' )" --args "/home/bgudenas/Proj/${PROJ}/Results/${BASE}_CountMat.csv" "/home/bgudenas/Proj/${PROJ}/Data/${BASE}/${BASE}_STAR_Logs.csv"
fi
