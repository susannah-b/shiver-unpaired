#!/usr/bin/env bash

################################################################################
# What do you have to type into the command line to make these commands execute?
# (If the binary file lives in a directory that is not included in your $PATH
# variable, you will need to include the path here.)
BlastDBcommand='makeblastdb'
smalt='smalt'
samtools='samtools'
trimmomatic="$HOME/Trimmomatic-0.33/trimmomatic-0.33.jar"
# Newer versions of fastaq have a space instead of an underscore between 
# 'fastaq' and the desired fastaq command, e.g. 'sequence_trim'. This means we
# will interpret space in this command as meaningful (i.e. we will execute
# $FastaqSequenceTrim not "$FastaqSequenceTrim"), so if you are including a path
# to the fastaq binary, it may not contain whitespace.
FastaqSequenceTrim='fastaq_sequence_trim'

# The trimmomatic manual explains at length the parameters controlling read
# trimming; the reader is referred to it for explanations of the following two
# variables and other options not used here:
IlluminaClipParams='2:10:7:1:true'
BaseQualityParams='MINLEN:50 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20'
NumThreadsTrimmomatic=1
smaltIndexOptions="-k 15 -s 3"
smaltMapOptions="-x -y 0.75 -j 0 -i 2000"
samtoolsReadFlags='-f 3 -F 4'
mpileupOptions='--min-BQ 5 --max-depth 1000000'

MinCov1=10
MinCov2=20

# Suffixes we'll append to the sample ID for output files.
# If you change the extension, you may well break something.
OutputRefSuffix='_ref.fasta'
MappedContaminantReadsSuffix='_ContaminantReads.bam'
BaseFreqsSuffix='_BaseFreqs.csv'
InsertSizeCountsSuffix='_InsertSizeCounts.csv'
CoordsDictSuffix='_coords.csv'
CleanedReads1Suffix='_clean_1.fastq' # .gz will be added when they're zipped
CleanedReads2Suffix='_clean_2.fastq' # .gz will be added when they're zipped
GlobalAlnSuffix='_ForGlobalAln.fasta'
################################################################################
# The names of temporary files we'll create in the working directory.
# If you change the extension, you may well break something.
TempRefAlignment='temp_RefAlignment.fasta'
GappyRefWithExtraSeq='temp_GappyRefWithExtraSeq.fasta'
FlattenedContigs='temp_FlattenedContigs.fasta'
AllContigsList='temp_AllContigsList.txt'
HIVContigsList='temp_HIVContigsList.txt'
ContaminantContigsList='temp_ContaminantContigsList.txt'
RefAndContaminantContigs='temp_RefAndContaminantContigs.fasta' # no whitespace!
BlastDB='temp_BlastDB' # no whitespace!
BadReadsBaseName='temp_ContaminantReads'
smaltIndex='temp_smaltRefIndex'
AllMappedContaminantReads='temp_ContaminantReads_IncUnmapped.sam'
RefFromAlignment='temp_RefFromAlignment.fasta'
AllSeqsInAln='temp_AllSeqsInAln.txt'
reads1asFasta='temp_reads1.fasta'
reads2asFasta='temp_reads2.fasta'
reads1blast1='temp_reads1_1.blast'
reads2blast1='temp_reads2_1.blast'
reads1blast2='temp_reads1_2.blast'
reads2blast2='temp_reads2_2.blast'
reads1sorted='temp_1_sorted.fastq'
reads2sorted='temp_2_sorted.fastq'
MapOutAsSam='temp_MapOut.sam'
MapOutConversion1='temp_MapOutStep1.bam'
MapOutConversion2='temp_MapOutStep2'
MapOutConversion3='temp_MapOutStep3.bam'
InsertSizes1='temp_InsertSizes.txt'
InsertSizes2='temp_InsertSizes2.txt'
PileupFile='temp_MapOut.pileup'
RefWithGaps='temp_RefWithGaps.fasta'
reads1trim1='temp_reads1trim1.fastq'
reads1trim2='temp_reads1trim2.fastq'
reads2trim1='temp_reads2trim1.fastq'
reads2trim2='temp_reads2trim2.fastq'
reads1trimmings='temp_trimmings1.fastq'
reads2trimmings='temp_trimmings2.fastq'
AlignmentForTesting='temp_test.fasta'
