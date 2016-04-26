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

# Shall we trim (remove adapters, primers and low quality bases)?
# Shall we clean (remove read pairs that look like contaminants)?
TrimReads=true
CleanReads=true

# The trimmomatic manual explains at length the parameters controlling read
# trimming; the reader is referred to it for explanations of the following two
# variables and other options not used here:
IlluminaClipParams='2:10:7:1:true'
BaseQualityParams='MINLEN:50 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20'
NumThreadsTrimmomatic=1

# Check the smalt documentation for a full explanation of mapping parameters,
# including those not used by default here.
# A summary of the index options used below:
# -k sets the word (kmer) length, -s the sampling step size (i.e. is every word
#  hashed, or every second word, or one word in every 3, ...), when a hash table
# is made using the reference.
# A summary of the mapping options used below:
# -x means a read and its mate are mapped independently (not constraining them
# to be close), -y sets the minimum fraction of identical nucleotides a read
# must have to its reference before it is considered mapped, -j is the minimum 
# insert size and -i the maximum insert size: outside of this range, the read
# pair is still mapped, but flagged as improperly paired.
smaltIndexOptions="-k 15 -s 3"
smaltMapOptions="-x -y 0.75 -j 0 -i 2000"

# For a description of the SAM format flags, look at
# https://samtools.github.io/hts-specs/SAMv1.pdf
# Note that the flags are combined in a bitwise manner, which is fiddly. See 
# http://onetipperday.blogspot.co.uk/2012/04/understand-flag-code-of-sam-format.html
# for a discussion. The flags used below mean unmapped reads are excluded (-F 4)
# and only properly aligned pairs are kept (-f 3).
samtoolsReadFlags='-f 3 -F 4'

# See http://www.htslib.org/doc/samtools.html for a description of samtools
# mpileup options. Those used below mean that the minimum of the base quality
# the 'BAQ' quantity (see http://samtools.sourceforge.net/mpileup.shtml for an
# explanation) must be at least 5, and only the first 1000000 reads mapped to
# each point will be considered (NB a limit must be provided; the default is
# 250).
mpileupOptions='--min-BQ 5 --max-depth 1000000'

# The minimum number of reads mapped to reference to call a base at that
# position: if there are fewer reads we call '?' instead of a base.
MinCov1=10
# The minimum number of reads to use upper case (to signal increased confidence)
# for the base called.
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
RawContigFile='temp_HIVcontigs_raw.fasta'
CutContigFile='temp_HIVcontigs_cut.fasta'
TempRawContigAlignment='temp_HIVcontigs_raw_wRefs_swap.fasta'
TempCutContigAlignment='temp_HIVcontigs_cut_wRefs_swap.fasta'
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
