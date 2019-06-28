#BSUB -P PDX
#BSUB -R "rusage[mem=10000]"
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -J "star[1-8]"
set -ue pipefail

RUN="YES"
FDIR=./Raw
NCUT=25
HCUT=9
THR=1

FASTQS=$( find ${FDIR}* -type f -name "*q.gz" -exec ls {} + )
SAMP=$( ls $FASTQS | rev | cut -c ${NCUT}- | rev  |sed 's!.*/!!' | sed 's/_$//' | sort | uniq -d | awk -v var="$LSB_JOBINDEX" 'FNR == var {print}' )
export PATH=/home/bgudenas/src/:$PATH

#1. Adapter trim
#2. Map reads mouse and human
#3. Disambiguate reads

GetReads.sh $FDIR $NCUT $HCUT $LSB_JOBINDEX $RUN

if [ "$RUN" = "YES" ]
	then
AdapterTrim.sh $THR
StarPDX.sh $THR

fi

