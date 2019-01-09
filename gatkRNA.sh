#BSUB -P Clingen
#BSUB -R "rusage[mem=30000]"
#BSUB -J "Clingen[1-9]"
#BSUB -n 1

RUN="YES"
BASE="Clingen"
NCUT=30

## USE $LSB_JOBINDEX to select sample from SHHs.txt
BAM=$( ls ./BAM/${BASE}/*.bam | awk -v ind="$LSB_JOBINDEX" 'FNR == ind {print}' )
SAMP=$( basename $BAM | rev | cut -c ${NCUT}- | rev )

RUN="YES"
PLATFORM="illumina"
MACHINE="Hiseq2500"

echo $TMPDIR

module purge
module load gatk/3.7
module load picard/2.9.4
module load annovar/043014
module load samtools

PICARD=/hpcf/apps/picard/install/2.9.4/picard.jar
GATK=/hpcf/apps/gatk/install/3.7/GenomeAnalysisTK.jar

## Call Variants using GATK RNA-seq protocol
# BAMs should be in "./BAM/"
### References
FASTA="/home/bgudenas/Annots/Human/Homo_Sapiens_hg38_refChromes_ensembl.fa"
REFSNP="/home/bgudenas/Annots/Human/VCF/ensembl_common_all_hg38p7.vcf"
INTERVAL="/home/bgudenas/Proj/SHHalpha/ELP1.intervals"
INTERVAL_LIST="/home/bgudenas/Proj/SHHalpha/RefChromes_hg38_99.interval_list"
HUMANDB="/home/bgudenas/Proj/SHHalpha/humandb/"

mkdir -p ./GATK/${BASE}/VCF/Annovar

    printf "\n"
    echo "SAMPLE##################"
    echo $SAMP
    

    echo "BAM---####"
    ls -lh $BAM

if [ "$RUN" = "YES" ]
	then

	if [ ! -f ./GATK/${BASE}/${SAMP}_MC.bam ]
	 	then
java -jar -XX:ParallelGCThreads=1 \
	$PICARD FixMateInformation \
	I=$BAM \
	O=./GATK/${BASE}/${SAMP}_MC.bam \
	ADD_MATE_CIGAR=true \
	SO=coordinate \
	VALIDATION_STRINGENCY=LENIENT \
	TMP_DIR=$TMPDIR
fi

#if [ ! -f ./GATK/${BASE}/${SAMP}_SO.bam ]
#	then
#java -jar -Djava.io.tmpdir=$TMPDIR -XX:ParallelGCThreads=1 \
#	$PICARD SortSam \
#	SO=coordinate \
#        I=./GATK/${BASE}/${SAMP}_MC.bam \
#        O=./GATK/${BASE}/${SAMP}_SO.bam 

#rm ./GATK/${BASE}/${SAMP}_MC.bam
#fi

if [ ! -f ./GATK/${BASE}/${SAMP}_RG.bam ]
                then

java -jar -Djava.io.tmpdir=$TMPDIR -XX:ParallelGCThreads=1 \
	$PICARD FilterSamReads \
	I=./GATK/${BASE}/${SAMP}_MC.bam \
	O=./GATK/${BASE}/${SAMP}_region.bam \
	INTERVAL_LIST=$INTERVAL_LIST \
        FILTER=includePairedIntervals \
	VALIDATION_STRINGENCY=LENIENT

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

			## Remove alt contigs from BAM to match 1000g SNP VCF
		samtools view -H ./GATK/${BASE}/${SAMP}_dup.bam | grep -v "SN:KI" | grep -v "SN:GL" | grep -v "MT" | samtools reheader - ./GATK/${BASE}/${SAMP}_dup.bam > ./GATK/${BASE}/${SAMP}_new_head.bam

## Reorder Bam
		java -jar -XX:ParallelGCThreads=1 \
		    $PICARD ReorderSam \
		    I=./GATK/${BASE}/${SAMP}_new_head.bam \
		    O=./GATK/${BASE}/${SAMP}_reordered.bam \
		    R=$FASTA \
		    CREATE_INDEX=TRUE

java -jar -Djava.io.tmpdir=$TMPDIR -XX:ParallelGCThreads=1 \
		$GATK -T SplitNCigarReads \
		-R $FASTA \
		-I ./GATK/${BASE}/${SAMP}_reordered.bam \
		-o ./GATK/${BASE}/${SAMP}_split.bam \
		-rf ReassignOneMappingQuality \
		-RMQF 255 -RMQT 60  \
		-U ALLOW_N_CIGAR_READS
	
		fi

		if [ ! -f ./GATK/${BASE}/${SAMP}_BQSR.bam ]
	 	then
java -jar -Djava.io.tmpdir=$TMPDIR -XX:ParallelGCThreads=1 \
		$GATK -T BaseRecalibrator \
			-R $FASTA \
			-knownSites $REFSNP \
			-I ./GATK/${BASE}/${SAMP}_split.bam \
			-o ./GATK/${BASE}/${SAMP}_recal.table

java -jar -Djava.io.tmpdir=$TMPDIR -XX:ParallelGCThreads=1 \
		$GATK -T PrintReads \
			-R $FASTA \
			-BQSR ./GATK/${BASE}/${SAMP}_recal.table \
			-I ./GATK/${BASE}/${SAMP}_split.bam \
			-o ./GATK/${BASE}/${SAMP}_BQSR.bam
	fi

if [ ! -f ./GATK/${BASE}/VCF/${SAMP}_HC.vcf ]
	 	then

java -jar -Djava.io.tmpdir=$TMPDIR -XX:ParallelGCThreads=1 \
		$GATK -T HaplotypeCaller \
		-R $FASTA \
		-D $REFSNP \
		-I ./GATK/${BASE}/${SAMP}_BQSR.bam \
		-dontUseSoftClippedBases \
		-stand_call_conf 25.0  \
		-o ./GATK/${BASE}/VCF/${SAMP}_HC.vcf

		fi

if [ ! -f ./GATK/${BASE}/VCF/${SAMP}_Filt.vcf ]
	 	then

java -jar -Djava.io.tmpdir=$TMPDIR -XX:ParallelGCThreads=1 \
		$GATK -T VariantFiltration \
		-R $FASTA \
		-V ./GATK/${BASE}/VCF/${SAMP}_HC.vcf \
		-window 35 \
		-cluster 3 \
		-filterName FS \
		-filter "FS > 40.0" \
		-filterName QD \
		-filter "QD < 2.0" \
		-o ./GATK/${BASE}/VCF/${SAMP}_Filt.vc

		fi

convert2annovar.pl ./GATK/${BASE}/VCF/${SAMP}_Filt.vc -format vcf4 > ./GATK/${BASE}/VCF/${SAMP}_annot 
table_annovar.pl ./GATK/${BASE}/VCF/${SAMP}_annot $HUMANDB -buildver hg38 \
-out ./GATK/${BASE}/Annovar/${SAMP} -remove -protocol ensGene,cytoBand,exac03,avsnp150,dbnsfp35c -operation g,r,f,f,f -nastring NA --csvout

fi

