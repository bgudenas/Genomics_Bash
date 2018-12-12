#BSUB -P STIE
#BSUB -R "rusage[mem=3000]"
#BSUB -n 16
#BSUB --o mapped.out --eo mapped.err

RUN="NO"
NCUT=12
BASE="SRA"
THR=16

## Stringtie - transcript quantification

### BAM files should be in ./Data/${BASE}/ with *XS.bam


### ENSEMBL annots
GTF="../../Annots/Mouse/mmGenes92.gtf"

module purge
printf "STRINGTIE ASSEMBLY ############ \n"
date
stringtie --version
module load R/3.4.0



#########################################

SAMPNUM=$(ls ./Data/${BASE}/*XS.bam  | wc -l)
printf "SAMPLES DETECTED = ${SAMPNUM} \n"

mkdir -p ./Stringtie/${BASE}

for i in $(ls ./Data/${BASE}/*XS.bam | rev | cut -c 7- | rev)
do	
    printf "\n"  
    SAMP=$(basename ${i})
    echo "SAMPLE##################"
    echo $SAMP
    
    BAM=${i}"XS.bam"
    echo "BAM---####"
    ls -lh $BAM
    OUT="./Stringtie/"${BASE}"/"${SAMP}".gtf"

if [ "$RUN" = "YES" ]
	then

  stringtie -p $THR --rf -G $GTF -o $OUT -l $SAMP $BAM
fi
done

if [ "$RUN" = "YES" ]
  then

printf "######################################\n"
printf "MERGE\n"
printf "######################################"

mkdir -p ./Stringtie/${BASE}/Assembled

MERGE="./Stringtie/"${BASE}"/Assembled/merge_list.txt"
ls ./Stringtie/${BASE}/*.gtf > $MERGE

cat $MERGE
OUTGTF="./Stringtie/"${BASE}"/Assembled/stringtie_merged.gtf"

stringtie --merge -p $THR -F 0.5  -m 100 -G $GTF -o $OUTGTF $MERGE 
awk '$3=="transcript"' $OUTGTF | wc -l

printf "######################################\n"
printf "BALLGOWN\n"
printf "######################################"

mkdir -p ./Stringtie/${BASE}/Ballgown

for i in $(ls ./Data/${BASE}/*XS.bam | rev | cut -c 7- | rev)
do  
    printf "\n"  
    SAMP=$(basename ${i})
    echo "SAMPLE##################"
    echo $SAMP
    
    BAM=${i}"XS.bam"

    GENE_ABS="./Stringtie/"${BASE}"/Ballgown/"${SAMP}"/gene_abundances.tsv"
    OUT="./Stringtie/"${BASE}"/Ballgown/"${SAMP}"/"${SAMP}".gtf"

  stringtie -p $THR -B -e --rf -G $OUTGTF -A $GENE_ABS -o $OUT $BAM
done


cd ./Stringtie/${BASE}/Ballgown/
Rscript /home/bgudenas/src/gene_expression_matrix.R
cd ../../../

fi

date
