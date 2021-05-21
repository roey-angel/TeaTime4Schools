#!/bin/bash

# Description: run FUNGuild/NAMEGuild on ITS OTU table
# Instructions: go to the Analysis folder and run:
# 02_funguild.sh DATAFOLDER


## TODO

logFile="./02_funguild.log"

HOMEFOLDER=$(pwd)
DB="fungi"

touch ${HOMEFOLDER}/${logFile}
echo -e "02_funguild.sh  \n" >${HOMEFOLDER}/${logFile}

if [ -z ${1+x} ] # input argument?
    then 
        DATAFOLDER="../../Data/ITS/MiSeq_431_Roey__Angel_24541_2_-126504387/noPrimers/DADA2_pseudo/"
        echo "DATAFOLDER not provided. Will use ../../Data/ITS/MiSeq_431_Roey__Angel_24541_2_-126504387/noPrimers/DADA2_pseudo/" >> ${HOMEFOLDER}/${logFile}
    else 
        DATAFOLDER=$1
        echo "DATAFOLDER is set to '$1'" >> ${HOMEFOLDER}/${logFile}
fi

paste ${DATAFOLDER}/DADA2.seqtab_nochim.tsv <(awk 'BEGIN {printf "taxonomy\n"} NR>1 {for(i=2;i<=8;++i) printf $i";" ; print ""}' ${DATAFOLDER}/DADA2.taxa_unite.tsv) > ${DATAFOLDER}/DADA2.seqtab_nochim_funguild.tsv

## Run funguild script
echo -e "Starting FUNGuild python script\n" >> ${HOMEFOLDER}/${logFile}
python /usr/local/bin/funguild -otu ${DATAFOLDER}/DADA2.seqtab_nochim_funguild.tsv -db $DB -m -u >> /${HOMEFOLDER}/${logFile} 2>&1
rm DADA2.seqtab_nochim_funguild.tsv

