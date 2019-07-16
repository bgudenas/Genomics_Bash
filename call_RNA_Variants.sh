#!/bin/bash

echo "#########"
echo "Call Variants"
echo 'Threads = ' $1
THR=$1
echo 'Read Length = ' $2
RL=$2
module purge
module load star/2.6.1c
module load sambamba

R1=$(ls $TMPDIR/*1.fastq )
ls $R1
R2=$(ls $TMPDIR/*2.fastq )
ls $R2

SAMP=$( basename $R1 | rev | cut -c 10- | rev )
echo $SAMP
mkdir -p $TMPDIR/BAM
mkdir -p ./GATK/BQSR
mkdir -p ./GATK/HC
mkdir -p ./GATK/Filt
mkdir -p ./GATK/Metrics

######### References
### GATK bundle
BUNDLE=/research/rgs01/project_space/northcgrp/Northcott_Bioinformatics/northcgrp/References/Human/GATK/hg38/
FASTA=$BUNDLE/Homo_sapiens_assembly38.fasta
REFSNP=$BUNDLE/dbsnp_146.hg38.vcf.gz
REFINDEL=$BUNDLE/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
HSINDEL=$BUNDLE/Homo_sapiens_assembly38.known_indels.vcf.gz
INTERVAL=$BUNDLE/wgs_calling_regions.hg38.interval_list

GENDIR=$BUNDLE/STAR_index
##################
GATK=/home/bgudenas/Software/gatk-4.1.1.0/gatk
module load java/1.8.0_181 

DUP=$TMPDIR/BAM/${SAMP}_DUP.bam
if [ ! -f $DUP ]
then

INBAM=${TMPDIR}/${SAMP}Aligned.out.bam
SORTBAM=${TMPDIR}/${SAMP}_Sorted.bam

STAR --runMode alignReads --genomeDir $GENDIR --sjdbOverhang $RL --readFilesIn $R1 $R2 \
--runThreadN $THR \
--outFileNamePrefix $TMPDIR/${SAMP} \
--outSAMtype BAM Unsorted \
--twopassMode Basic \
--limitSjdbInsertNsj 2000000 \
--limitOutSJcollapsed 2000000 \
--outTmpDir $TMPDIR/tmp${SAMP}

sambamba sort -t $THR -m 38G --tmpdir $TMPDIR \
-o $SORTBAM $INBAM

$GATK MarkDuplicates \
-I $SORTBAM \
-O $DUP \
-M "./GATK/Metrics/${SAMP}_metric.txt" \
--CREATE_INDEX TRUE 1> /dev/null
fi

SPLIT=$TMPDIR/BAM/${SAMP}_SPLIT.bam
if [ ! -f $SPLIT ]
then
$GATK AddOrReplaceReadGroups \
-I  $DUP \
-O $TMPDIR/${SAMP}_RG.bam \
-LB $SAMP \
-PL "illumina" \
-PU "Hiseq2500" \
--CREATE_INDEX TRUE \
-SM $SAMP

#samtools index $TMPDIR/${SAMP}_RG.bam

$GATK SplitNCigarReads \
-R $FASTA \
-I $TMPDIR/${SAMP}_RG.bam \
-L $INTERVAL \
--tmp-dir $TMPDIR \
-O $SPLIT
fi

FIN=$TMPDIR/BAM/${SAMP}_GATK.bam
if [ ! -f $FIN ]
then

$GATK BaseRecalibrator \
-R $FASTA \
--known-sites $REFSNP \
--known-sites $REFINDEL \
--known-sites $HSINDEL \
-L $INTERVAL \
-I $SPLIT \
-O "./GATK/BQSR/${SAMP}-recal-table.txt" 1> /dev/null

$GATK ApplyBQSR \
-R $FASTA \
-I $SPLIT \
-L $INTERVAL \
--bqsr-recal-file "./GATK/BQSR/${SAMP}-recal-table.txt" \
-O $FIN 1> /dev/null

fi

### Variant calling
VCF_HC=./GATK/HC/${SAMP}.vcf.gz
if [ ! -f $VCF_HC ];
then

$GATK HaplotypeCaller  \
-R $FASTA \
-I $FIN \
--bam-output ./GATK/HC/${SAMP}.bam \
--create-output-bam-index TRUE \
--dbsnp $REFSNP \
--dont-use-soft-clipped-bases \
-stand-call-conf 20.0 \
-L $INTERVAL \
--tmp-dir $TMPDIR \
-O $VCF_HC 
fi

MERGED=./GATK/Filt/${SAMP}_Merged.vcf.gz
if [ ! -f $MERGED ];
then

$GATK SelectVariants \
-V $VCF_HC \
-select-type SNP \
-O $TMPDIR/${SAMP}_SNP.vcf.gz

$GATK SelectVariants \
-V $VCF_HC \
-select-type INDEL \
-O $TMPDIR/${SAMP}_INDEL.vcf.gz

$GATK VariantFiltration \
-V $TMPDIR/${SAMP}_SNP.vcf.gz \
--cluster-size 3 \
--cluster-window-size 35 \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "FS > 30.0" --filter-name "FS30" \
-O $TMPDIR/${SAMP}_SNP_Filt.vcf.gz

$GATK VariantFiltration \
-V $TMPDIR/${SAMP}_INDEL.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "FS > 30.0" --filter-name "FS40" \
-O $TMPDIR/${SAMP}_INDEL_Filt.vcf.gz

$GATK MergeVcfs \
-I $TMPDIR/${SAMP}_SNP_Filt.vcf.gz \
-I $TMPDIR/${SAMP}_INDEL_Filt.vcf.gz \
-O $MERGED
fi

TRIM=./GATK/Filt/${SAMP}_Trim.vcf.gz
$GATK SelectVariants \
-V $MERGED \
-O $TRIM \
--exclude-filtered

ANNOVAR=/home/bgudenas/Software/annovar/
HUMANDB=/home/bgudenas/Software/annovar/humandb
mkdir ./Annovar

${ANNOVAR}/table_annovar.pl $TRIM /home/bgudenas/Software/annovar/humandb \
-buildver hg38 -out "Annovar/${SAMP}" \
-remove -protocol ensGene,refGene,cytoBand,exac03,dbnsfp35c -operation g,g,r,f,f -polish -nastring . --vcfinput

## ADD vcf label to annovar txt output to label normal-tumor format field
## also remove intronic variants
header=$( cat ./Annovar/${SAMP}.hg38_multianno.vcf  | grep "CHROM" | cut -c 2- )
module load R/3.4.0
Rscript /home/bgudenas/src/addlab.R ./Annovar/${SAMP}.hg38_multianno.txt $header

val=$( wc -l ./Annovar/${SAMP}.hg38_multianno.txt )
LINES=$( echo $val | grep -ohw "[0-9]*\+[[:space:]]\+" )
if [ $LINES -eq 1 ];
then
echo "NO VARIANTS"
  rm ./Annovar/${SAMP}.*i*
  fi
