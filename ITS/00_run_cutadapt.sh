#!/bin/bash

# Description: run cutadapt to remove primer sequences
# Instructions: go to the Analysis folder and run:
# 00_run_cutadapt.sh $DATAFOLDER $FWD $REV 
# where $FWD and $REV are the forward and reverse primer sequences
# Note: multicore support (-j) only works with python 3

## TODO

nCore=0 # use all cores
minReadLen=50 # minimal read length (don't set to 0!)
logFile="./00_cutadapt.log"

HOMEFOLDER=$(pwd)
SCRIPTS=${HOMEFOLDER}

if [ -z ${1+x} ] # input argument?
    then 
        logFile="cutadapt.log"
        DATAFOLDER="../../Data/ITS/MiSeq_431_Roey__Angel_24541_2_-126504387"
        #FWD="ACACTGACGACATGGTTCTACACTTGGTCATTTAGAGGAAGTAA" # CS tag was not removed
        FWD="CTTGGTCATTTAGAGGAAGTAA" # CS tag was removed
        REV="GCTGCGTTCTTCATCGATGC" # CS tag was removed
        #REV="TACGGTAGCAGAGACTTGGTCTGCTGCGTTCTTCATCGATGC" # CS tag was not removed
        touch ${HOMEFOLDER}/${logFile}
        echo -e "00_run_cutadapt.sh  \n" >${HOMEFOLDER}/${logFile}        
        echo "Primers not provided. Will use ITS1f - CTTGGTCATTTAGAGGAAGTAA, ITS2 - GCTGCGTTCTTCATCGATGC" >> ${HOMEFOLDER}/${logFile}
    else 
        DATAFOLDER=$1
        FWD=$2
        REV=$3
        logFile=$4
        touch ${HOMEFOLDER}/${logFile}
        echo -e "00_run_cutadapt.sh  \n" >${HOMEFOLDER}/${logFile}
        echo "FWD is set to '$2' and REV is set to '$3'" >> ${HOMEFOLDER}/${logFile}
fi

## Merge files of identical samples from repeated runs (if any). Otherwise just copy all .gz files to one desitnation
dest=data_files
mkdir -p ${DATAFOLDER}/$dest
find `ls -d ${DATAFOLDER}/*/ | grep -v "$dest"` -name "$dest" -prune -o -name '*.fastq.gz' | while read file
do    base=$(basename "$file")
      if [ -s "${DATAFOLDER}/$dest/$base" ]
      then cat "$file" # potentially manipulate the first file (e.g. sed 1d <"$file")
      else cat "$file"
      fi >>"${DATAFOLDER}/$dest/$base"
done

## Run cutadapt R script
echo -e "Starting R script \n" >> ${HOMEFOLDER}/${logFile}
/usr/bin/time -v Rscript --vanilla ${SCRIPTS}/00_cutadapt.R ${DATAFOLDER} ${minReadLen}  ${FWD} ${REV} >> /${HOMEFOLDER}/${logFile} 2>&1
rm -rf ${DATAFOLDER}/$dest
