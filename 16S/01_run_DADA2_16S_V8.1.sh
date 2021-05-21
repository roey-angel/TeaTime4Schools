#!/bin/bash
#

# Description: run DADA2 pipeline on a MiSeq run data
# Instructions: go to the Analysis folder and run:
# 01_run_DADA2_16S_V8.1.sh $POOLING $DATAFOLDER $RESOURCES

HOMEFOLDER=$(pwd)
POOLING=${1:-pseudo} # pseudo, TRUE, FALSE
DATAFOLDER=${2:-"./tmp_data"}
RESOURCES=${3:-"~/Resources/DADA2/"} # RESOURCES should be referred to from $DATAFOLDER/DADA2_pseudo!
SCRIPTS=${HOMEFOLDER}
GENE=16S
VERSION=V8.1

# Collect all samples
if [ -d $DATAFOLDER ]
then
 rm -rf $DATAFOLDER
fi
mkdir $DATAFOLDER
cp ../../1st_workshop/Data/16S/noPrimers/*.* $DATAFOLDER
cp ../../2nd_workshop/Data/16S/noPrimers/*.* $DATAFOLDER
cp ../../3rd_workshop/Data/16S/noPrimers/*.* $DATAFOLDER
cp ../../4th_workshop/Data/16S/noPrimers/*.* $DATAFOLDER

cd $DATAFOLDER
rm Eva28_Gerrit0_1_S237_L001_R1_001_noPrimers.fastq.gz Eva28_Gerrit0_1_S237_L001_R2_001_noPrimers.fastq.gz Eva31_Gerrit0_22_S240_L001_R1_001_noPrimers.fastq.gz Eva31_Gerrit0_22_S240_L001_R2_001_noPrimers.fastq.gz Eva34_GerritNTC_S243_L001_R1_001_noPrimers.fastq.gz Eva34_GerritNTC_S243_L001_R2_001_noPrimers.fastq.gz
cd ../

if [ -d $DATAFOLDER/DADA2_pseudo ]
then
 rm -rf $DATAFOLDER/DADA2_pseudo
fi
mkdir $DATAFOLDER/DADA2_pseudo


LOG=01_run_DADA2_${GENE}_${POOLING}_${VERSION}.log
touch ${HOMEFOLDER}/${LOG}

echo -e "01_run_DADA2_${POOLING}_${GENE}_${VERSION}.sh  \n" >${HOMEFOLDER}/${LOG}
echo -e "Will process reads from the following folder:"  >> ${HOMEFOLDER}/${LOG}
echo -e $DATAFOLDER >> ${HOMEFOLDER}/${LOG}

echo -e "Starting R script \n" >> ${HOMEFOLDER}/${LOG}
cd $DATAFOLDER/DADA2_${POOLING}
/usr/bin/time -v Rscript --vanilla ${SCRIPTS}/01_DADA2_16S_${VERSION}.R $POOLING ../ $RESOURCES >> ${HOMEFOLDER}/${LOG} 2>&1
cd ${HOMEFOLDER}
mv $DATAFOLDER/DADA2_${POOLING} ${HOMEFOLDER}
rm -rf ${DATAFOLDER}
