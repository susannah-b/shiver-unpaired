#!/usr/bin/env bash

################################################################################
BlastDBcommand='makeblastdb'
smalt='smalt'
samtools='samtools'

smaltIndexOptions="-k 15 -s 3"
smaltMapOptions="-x -y 0.75 -j 0 -i 2000"
samtoolsReadFlags='-f 3 -F 4'
mpileupOptions='--min-BQ 5 --max-depth 1000000'

MinCoverage1=10
MinCoverage2=20

################################################################################
# The names of some temporary files we'll create in the working directory
TempRefAlignment='temp_RefAlignment.fasta'
GappyRefWithExtraSeq='temp_GappyRefWithExtraSeq.fasta'
FlattenedContigs='temp_FlattenedContigs.fasta'
AllContigsList='temp_AllContigsList.txt'
HIVContigsList='temp_HIVContigsList.txt'
ContaminantContigsList='temp_ContaminantContigsList.txt'
RefAndContaminantContigs='temp_RefAndContaminantContigs.fasta'
BlastDB='temp_BlastDB'
BadReadsBaseName='temp_ContaminantReads'
smaltIndex='temp_smaltRefIndex'
AllMappedContaminantReads='ContaminantReads_IncUnmapped.sam'

