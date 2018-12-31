#BSUB -P Stringtie
#BSUB -R "rusage[mem=10000]"
#BSUB -n 5
#BSUB -J "Stringstar[1-170]"
#BSUB -R "span[hosts=1]"

#set -ueo pipefail 

RUN="YES"
MERGE="NO"
TWOPASS="NO"
THR=5

module load samtools/1.9

## USE $LSB_JOBINDEX to select sample from SHHs.txt
SAMP=$( ls ./Raw | awk -v ind="$LSB_JOBINDEX" 'FNR == ind {print}' )

GTF="/home/bgudenas/Annots/Human/Homo_sapiens.GRCh38.93.gtf"

date
stringtie --version

    printf "\n"  
    echo "SAMPLE##################"
    echo $SAMP

BAM="/scratch_space/bgudenas/Landscape/"${SAMP}"XS.bam"

if [ ! -f $BAM ];
    then 

    BAM="./Data/"${i}"Aligned.sortedByCoord.out.bam"
    echo "BAM---####"
    ls -lh $BAM

## Add XS strand to BAM
samtools view --threads $THR -h $BAM | awk -v strType=2 -f /home/bgudenas/src/tagXSstrandedData.awk | samtools view --threads $THR -bS - > "/scratch_space/bgudenas/Landscape/"${SAMP}"XS.bam"

fi

BAM="/scratch_space/bgudenas/Landscape/"${SAMP}"XS.bam"

OUTSTRING="./Stringtie/"${SAMP}".gtf"

if [ "$RUN" = "YES" ] && [ ! -f $OUTSTRING ];
    then

  stringtie -p $THR -j 2 --rf -G $GTF -o $OUTSTRING -l $SAMP $BAM
fi


OUTGTF="./Stringtie/Assembled/stringtie_merged.gtf"

if [ "$MERGE" = "YES" ] && [ ! -f $OUTGTF ]
  then

    printf "######################################\n"
    printf "MERGE\n"
    printf "######################################"

    mkdir -p ./Stringtie/Assembled

    MERGE="./Stringtie/Assembled/merge_list.txt"
    ls ./Stringtie/*.gtf > $MERGE

    cat $MERGE

    stringtie --merge -p $THR -G $GTF -o $OUTGTF $MERGE 
    awk '$3=="transcript"' $OUTGTF | wc -l

fi

GENE_ABS="./Stringtie/Ballgown/"${SAMP}"/gene_abundances.tsv"

if [ "$TWOPASS" = "YES" ] && [ ! -f $GENE_ABS ]
then

    OUT="./Stringtie/Ballgown/"${SAMP}"/"${SAMP}".gtf"

    stringtie -p $THR -B -e --rf -j 2 -G $OUTGTF -A $GENE_ABS -o $OUT $BAM

fi

  
