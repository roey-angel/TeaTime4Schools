#!/bin/bash

# Description: a wrapper to run 00_cutadapt.sh on several datasets
# Instructions: go to the analysis folder and run  00_run_cutadapt.sh

## TODO

bash 00_run_cutadapt.sh ../../1st_workshop/Data/ITS/ CTTGGTCATTTAGAGGAAGTAA GCTGCGTTCTTCATCGATGC cutadapt_1st.log &
bash 00_run_cutadapt.sh ../../2nd_workshop/Data/ITS/ CTTGGTCATTTAGAGGAAGTAA GCTGCGTTCTTCATCGATGC cutadapt_2nd.log &
bash 00_run_cutadapt.sh ../../3rd_workshop/Data/ITS/ CTTGGTCATTTAGAGGAAGTAA GCTGCGTTCTTCATCGATGC cutadapt_3rd.log &
bash 00_run_cutadapt.sh ../../4th_workshop/Data/ITS/ CTTGGTCATTTAGAGGAAGTAA GCTGCGTTCTTCATCGATGC cutadapt_4th.log &
wait
