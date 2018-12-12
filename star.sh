#BSUB -P star_north
#BSUB -R "rusage[mem=3000]"
#BSUB -n 16
#BSUB --o mapped.out --eo mapped.err

RUN="YES"
NCUT=21
BASE="northcgrp_125169_totalstranded"
THR=16

## Map reads using STAR multi-sample protocol
## first map using basic params
## Store SJ junctions
## Second pass using all sample junctions

### fastq files should be in ./Raw/${BASE}/*/ with sample *.fastq.gz inside

### ENSEMBL annots
GENDIR="../../Annots/Mouse/mm10_star99"
GTF="../../Annots/Mouse/mmGenes92.gtf"

module purge
module load samtools/1.7
printf "STAR RNA-seq mapping############# \n"
date
STAR --version

rm outStar.log
touch outStar.log
#########################################

SAMPNUM=$(ls ./Raw/${BASE}*/*/*.fastq.gz | rev | cut -c ${NCUT}- | rev  |  sed 's!.*/!!' | sort | uniq -d | wc -l)
printf "SAMPLES DETECTED = ${SAMPNUM} \n"
printf "FASTQ'S DETECTED = $(ls -l ./Raw/${BASE}*/*/* | wc -l) \n"

RLENGTH=$(gunzip -c "$(ls ./Raw/${BASE}*/*/*.gz | head -n1)" | head -n 1000 | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c)
printf "READ LENGTH(n=250) = ${RLENGTH} \n\n"


for i in $(ls ./Raw/${BASE}*/*/*.fastq.gz | rev | cut -c ${NCUT}- | rev  |sed 's!.*/!!' | sort | uniq -d)
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
    MR1="./Data/Read1_${SAMP}.fastq.gz"
    MR2="./Data/Read2_${SAMP}.fastq.gz"
    
if [ "$RUN" = "YES" ]
	then
    cat $R1 > $MR1
    cat $R2 > $MR2
    
    ls -lh $MR1
    ls -lh $MR2

    echo "STAR"
    STAR --runMode alignReads --genomeDir $GENDIR --readFilesIn $MR1 $MR2\
    --runThreadN $THR --readFilesCommand gunzip -c \
    --outFileNamePrefix Data/${SAMP} \
    --outFilterType BySJout \
    --alignIntronMin 20 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax  0.05 \
    --outFilterMultimapNmax 20 \
    --peOverlapNbasesMin 10 \
    >> outStar.log 2>&1

    OUTBAM="./Data/"${SAMP}"Aligned.out.sam"
    rm $OUTBAM
fi
done

##########################################################################################################
## Second pass using all indexes --sjdbFileChrStartEnd
## Move junction files to prevent overwriting
if [ "$RUN" = "YES" ]
	then
mkdir ./Data/${BASE}

mkdir ./Data/${BASE}/SJ
mv ./Data/*SJ.out* ./Data/${BASE}/SJ

awk -f /home/bgudenas/src/sjCollapseSamples.awk ./Data/${BASE}/SJ/*SJ.out.tab | sort -k1,1V -k2,2n -k3,3n > ./Data/${BASE}/SJ/SJ.all
grep -E -v "KI|GL|MT" ./Data/${BASE}/SJ/SJ.all > ./Data/${BASE}/SJ/SJ.pass1
awk '{ if ($10 >= 2) { print } }' ./Data/${BASE}/SJ/SJ.pass1 > ./Data/${BASE}/SJ/SJ.Filtered

mkdir -p ./BAM/${BASE}
	fi

for i in $(ls ./Raw/${BASE}*/*/*.fastq.gz | rev | cut -c ${NCUT}- | rev  |sed 's!.*/!!' | sort | uniq -d)
do
    SAMP=$(basename ${i})
    echo "ROUND 2 -- ${SAMP}"
    MR1="./Data/Read1_"${SAMP}".fastq.gz"
    MR2="./Data/Read2_"${SAMP}".fastq.gz"
if [ "$RUN" = "YES" ]
	then
STAR --runMode alignReads --genomeDir $GENDIR --readFilesIn $MR1 $MR2 \
    --runThreadN $THR --outBAMsortingThreadN $THR --readFilesCommand gunzip -c \
    --outFileNamePrefix Data/${BASE}/${SAMP} \
    --outSAMtype BAM SortedByCoordinate --outSAMattrIHstart 0 \
    --quantMode GeneCounts  \
    --sjdbFileChrStartEnd ./Data/${BASE}/SJ/SJ.Filtered ./Data/${BASE}/SJ/${SAMP}SJ.out.tab \
    --outReadsUnmapped Fastx \
    --peOverlapNbasesMin 10 \
    --limitSjdbInsertNsj 1200000 \
    --outFilterType BySJout \
    --alignIntronMin 20 \
    --alignSJDBoverhangMin 1 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax  0.05 \
    --outFilterMultimapNmax 20 \
    --limitBAMsortRAM 50000000000 \
    >> outStar.log 2>&1
    
   #Add XS attribute
    OUTBAM="./Data/"${BASE}"/"${SAMP}"Aligned.sortedByCoord.out.bam"
   ls $OUTBAM
   samtools view --threads $THR -h $OUTBAM | awk -v strType=2 -f /home/bgudenas/src/tagXSstrandedData.awk | samtools view --threads $THR -bS - > "./Data/"${BASE}"/"${SAMP}"XS.bam"
   rm $OUTBAM
   rm $MR1
   rm $MR2
   rm -rf ./Data/${SAMP}*
   mkdir -p ./BAM/${BASE}
   mv ./Data/${BASE}/*XS.bam ./BAM/${BASE}/
   mv ./Data/${BASE}/*Unmapped* ./BAM/${BASE}/
	fi
done

cat ./Data/${BASE}/*Log.final* | grep -E 'Number of input reads|Uniquely mapped reads %'
date

module load R/3.4.0
cd ./Data/${BASE}
Rscript /home/bgudenas/src/countMat.R 
cd ../..
