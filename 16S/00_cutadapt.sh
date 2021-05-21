#!/bin/bash

# Description: run cutadapt to remove primer sequences
# Instructions: go to the Analysis folder and run:
# 00_cutadapt.sh $DATAFOLDER $FWD $REV
# where $FWD and $REV are the forward and reverse primer sequences
# Note: multicore support (-j) only works with python 3

## TODO

nCore=0 # use all cores
minReadLen=100 # minimal read length (do not set to 0!)
#logFile="./cutadapt.log"

HOMEFOLDER=$(pwd)

if [ -z ${1+x} ] # input argument?
    then 
        logFile="_00cutadapt.log"
        DATAFOLDER="../../1st_workshop/Data/16S/"
        #FWD="ACACTGACGACATGGTTCTACAGTGYCAGCMGCCGCGGTAA" # CS tag was not removed
        FWD="GTGYCAGCMGCCGCGGTAA" # CS tag was removed
        REV="GGACTACNVGGGTWTCTAAT" # CS tag was removed
        #REV="TACGGTAGCAGAGACTTGGTCTGGACTACNVGGGTWTCTAAT" # CS tag was not removed
        touch ${HOMEFOLDER}/${logFile}
        echo -e "00_run_cutadapt.sh  \n" >${HOMEFOLDER}/${logFile}
        echo "Primers not provided. Will use 515F_mod - GTGYCAGCMGCCGCGGTAA, 806r - GGACTACNVGGGTWTCTAAT" >> ${HOMEFOLDER}/${logFile}
    else 
        DATAFOLDER=$1
        FWD=$2
        REV=$3
        logFile=$4
        touch ${HOMEFOLDER}/${logFile}
        echo -e "00_run_cutadapt.sh  \n" >${HOMEFOLDER}/${logFile}
        echo "FWD is set to '$2' and REV is set to '$3'" >> ${HOMEFOLDER}/${logFile}
fi

# make dir
if [ -d ${DATAFOLDER}/noPrimers ]                                                                                                                                                                                                                                                                                                                                                         
    then
 rm -rf ${DATAFOLDER}/noPrimers
fi
mkdir ${DATAFOLDER}/noPrimers

## Merge files of identical samples from repeated runs (if any)
# This assumes that the read files are in folders, 1 per sample
dest=data_files
mkdir -p ${DATAFOLDER}/$dest
find `ls -d ${DATAFOLDER}/*/ | grep -v "$dest"` -name "$dest" -prune -o -name '*.fastq.gz' | while read file
do    base=$(basename "$file")
      if [ -s "${DATAFOLDER}/$dest/$base" ]
      then cat "$file" # potentially manipulate the first file (e.g. sed 1d <"$file")
      else cat "$file"
      fi >>"${DATAFOLDER}/$dest/$base"
done

## Run cutadapt
#for R1file in *_R1_*; do
#    cutadapt -g $FWD -o ./noPrimer/${R1file%.fastq.gz}_noPrimer.fastq.gz ${R1file} --discard-untrimmed >> $logFile;
#done

#for R2file in *_R2_*; do
#    cutadapt -g $REV -o ./noPrimer/${R2file%.fastq.gz}_noPrimer.fastq.gz ${R2file} --discard-untrimmed >> .$logFile;
#done

for R1file in ${DATAFOLDER}/${dest}/*_R1_*; do
    R2file=$(echo $R1file | sed 's/_R1_/_R2_/g');
    ~/.local/bin/cutadapt -j $nCore -g $FWD -G $REV --minimum-length $minReadLen --discard-untrimmed -o ${DATAFOLDER}/noPrimers/`basename ${R1file%.fastq.gz}_noPrimers.fastq.gz` -p ${DATAFOLDER}/noPrimers/`basename ${R2file%.fastq.gz}_noPrimers.fastq.gz` ${R1file} ${R2file} >> ${HOMEFOLDER}/${logFile};
done
rm -rf ${DATAFOLDER}/$dest
