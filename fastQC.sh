#BSUB -P fQC
#BSUB -R "rusage[mem=2000]"
#BSUB -n 8
#BSUB --o fastqc.out --eo fastqc.err

BASE="SRA"

module purge
module load fastqc/0.11.5


for dir in ./Raw/*/
do
    dir=${dir%*/}
   dir=${dir##*/}

    if [[ $dir == ${BASE}* ]];
	then
	echo $dir
	mkdir -p ./QC/${dir}
	ls ./Raw/${dir}/*/*.fastq.gz
fastqc --threads 8 ./Raw/${dir}/*/*.fastq.gz -o ./QC/${dir}
	fi
done
