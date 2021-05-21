#!/bin/bash

# Description: align DADA2 output sequences using PyNast and build a FastTree tree
# Instructions: after funning DADA2 pipeline, run from the data folder:
# 05_calc_FastTree.sh $full_path_to_rep_seqs.fa 
# Where $rep_seqs.fa is the full path to a fasta file containing OTU representatives or ESVs 

## TODO


HOMEFOLDER=`pwd`
DBs="/home/angel/Resources/DADA2/"

# Start
logFile="05_calc_FastTree_pseudo.log"

if [ -z ${1+x} ]
    then 
        echo "No fasta file provided. Exiting" > ${HOMEFOLDER}/${logFile}
        exit 1 >> ${HOMEFOLDER}/${logFile}
    else 
        inputPath=$1
fi
DATAFOLDER="${inputPath%/*}"
INPUTSEQS=$(basename $inputPath)

echo -e "05_calc_FastTree.sh ${INPUTSEQS}" > ${HOMEFOLDER}/${logFile}
echo "Will process DADA2 sequences from the following file:" >> ${HOMEFOLDER}/${logFile}
echo $(ls $INPUTSEQS) >> ${HOMEFOLDER}/${logFile}

if [ -d ${DATAFOLDER}/FastTree ]
then
 rm -rf ${DATAFOLDER}/FastTree
fi
mkdir ${DATAFOLDER}/FastTree

cp ${INPUTSEQS} ${DATAFOLDER}/FastTree
cd ${DATAFOLDER}/FastTree
/usr/bin/time -v mothur "# align.seqs(candidate=${INPUTSEQS}, template=${DBs}/silva.seed_v132.align, flip=t, processors=16)" &>> ${HOMEFOLDER}/${logFile}
#mv ../${INPUTSEQS%.fa*}.align ./
/usr/bin/time -v FastTree -gtr -nt < ${INPUTSEQS%.fa*}.align > ${INPUTSEQS%.fa*}.tree >> ${HOMEFOLDER}/${logFile}

