#BSUB -P BWA
#BSUB -R "rusage[mem=3000]"
#BSUB -n 16
#BSUB -J "BWA_Pineo[1-10]"

#EX: LSB_JOBINDEX=2

RUN="YES"
PROCESS="NO"
NCUT=21
BASE="gajjagrp_144268"
THR=16

##NCUT removes trailing characters from FASTQ string to extract unique sample ID
## Map reads using BWA per read group and process for GATK4

module purge
module load bwa/0.7.4
module load samtools/1.7
module load bbmap
ADAPTERS=/hpcf/apps/bbmap/install/37.28/resources/adapters.fa

### fastq files should be in ./Raw/${BASE}/*/ with sample *.fastq.gz inside

### GATK bundle
FASTA=/research/rgs01/project_space/northcgrp/Northcott_Bioinformatics/northcgrp/References/Human/GATK/hg38/Homo_sapiens_assembly38.fasta
REFSNP=/research/rgs01/project_space/northcgrp/Northcott_Bioinformatics/northcgrp/References/Human/GATK/hg38/dbsnp_146.hg38.vcf
REFINDEL=/research/rgs01/project_space/northcgrp/Northcott_Bioinformatics/northcgrp/References/Human/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf


SAMPNUM=$(ls ./Raw/${BASE}*/*/*.fastq.gz | rev | cut -c ${NCUT}- | rev  |  sed 's!.*/!!' | sort | uniq -d | wc -l)

SAMP=$(ls ./Raw/${BASE}*/*/*.fastq.gz | rev | cut -c ${NCUT}- | rev  |sed 's!.*/!!' | sort | uniq -d | awk -v var="$LSB_JOBINDEX" 'FNR == var {print}' )


    printf "\n"  
    echo "SAMPLE##################"
    echo $SAMP
    R1="./Raw/${BASE}*/*/"${SAMP}"*R1*.fastq.gz"
    R2="./Raw/${BASE}*/*/"${SAMP}"*R2*.fastq.gz"
    echo "READ1####"
    ls -lh $R1
    echo "READ2####"
    ls -lh $R2

 ### Extract lane IDS Ex: L003
LANES=$(ls $R1 | grep -oh "L[0-9]*")
for i in $LANES
    do
        R1="./Raw/${BASE}*/*/"${SAMP}"*"${i}"_R1*.fastq.gz"
        R2="./Raw/${BASE}*/*/"${SAMP}"*"${i}"_R2*.fastq.gz"

###################################
    MR1="${TMPDIR}/Read1_${SAMP}.fastq"
    MR2="${TMPDIR}/Read2_${SAMP}.fastq"

if [ "$RUN" = "YES" ]
    then
    cat $R1 | gunzip > $MR1
    cat $R2 | gunzip > $MR2
    
    echo "LANE = "$i    
    ls -lh $MR1
    ls -lh $MR2

TR1="${TMPDIR}/Read1_Trim_${SAMP}.fastq"
TR2="${TMPDIR}/Read2_Trim_${SAMP}.fastq"

     bbduk.sh -Xmx32g threads=$THR \
        in1=$MR1 in2=$MR2 \
        out1=$TR1 out2=$TR2 \
        ref=$ADAPTERS ktrim=r ordered \
        k=22 mink=11 hdist=1 tbo tpe \
        rcomp=f \
        ow=t \
        minlen=30 \
        trimq=10 \
        qtrim=rl

bwa mem -aM -t $THR -v 2 \
-R "@RG\tID:${SAMP}${i}\tSM:${SAMP}\tPL:ILLUMINA\tLB:LIB" \
$FASTA \
$TR1 $TR2 | samtools sort -@ $THR - > $TMPDIR/${SAMP}_${i}.bam 

    fi
done

## MErge read groups bams
MERGED=/scratch_space/bgudenas/Pineo/${SAMP}merged.bam
ls -lh $TMPDIR/${SAMP}*
samtools merge -out $MERGED $TMPDIR/${SAMP}*bam


if [ "$PROCESS" == "YES" ]
    then

module purge
module load samtools/1.7
module load module load gatk/4.0.2.1

gatk MarkDuplicates \
        -I "$TMPDIR/${SAMP}merged.bam" \
        -O "/scratch_space/bgudenas/Pineo/GATK/${SAMP}dup.bam" \
        -M "./GATK/${SAMP}_metric.txt"

        gatk BaseRecalibrator \
                -R $FASTA \
                --known-sites $REFSNP \
                --known-sites $REFINDEL \
                -I "/scratch_space/bgudenas/Pineo/GATK/${SAMP}dup.bam" \
                -O "./GATK/${SAMP}-recal-table.txt"
            
            gatk ApplyBQSR \
                -R $FASTA \
                -I "./GATK/${SAMP}_dup.bam" \
                --bqsr-recal-file "./GATK/${SAMP}-recal-table.txt" \
                -O "/scratch_space/bgudenas/Pineo/GATK/${SAMP}GATK.bam"
    fi