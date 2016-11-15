#!/usr/bin/env bash

################################################################################
# Two options only needed for the fully automatic version: the maximum allowed
# percentage of gaps inside contigs when aligned to their closest reference (too
# much gap content indicates misalignment, rather than deletions), and the
# minimum fraction of a contig's length that blasts to HIV.
MaxContigGappiness=0.05
MinContigHitFrac=0.9

# What do you have to type into the command line to make these commands execute?
# (If the binary file lives in a directory that is not included in your $PATH
# variable, you will need to include the path here.)
BlastDBcommand='makeblastdb'
smalt='smalt'
samtools='samtools'
mafft='mafft'
# Newer versions of fastaq have a space instead of an underscore between 
# 'fastaq' and the desired fastaq command, e.g. 'sequence_trim'. This means we
# will interpret space in this command as meaningful (i.e. we will execute
# $FastaqSequenceTrim not "$FastaqSequenceTrim"), so if you are including a path
# to the fastaq binary, it may not contain whitespace.
FastaqSequenceTrim='fastaq sequence_trim'
# If you've downloaded the trimmomatic executable file (ending in .jar), to run
# it you probably need to type something like this:
# java -jar path/to/where/it/lives/trimmomatic-0.36.jar
# If someone else installed it for you (e.g. on MRC CLIMB) there may be an alias
# which means you just type 'trimmomatic' to run it:
trimmomatic="trimmomatic"

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
smaltMapOptions="-x -y 0.7 -j 0 -i 2000"

# After mapping, the choice of what kinds of reads should be kept is specified
# with SAM format flags, whose documentation is here:
# https://samtools.github.io/hts-specs/SAMv1.pdf
# Flags are combined in a bitwise manner, which is fiddly. This page
# https://broadinstitute.github.io/picard/explain-flags.html
# gives a more user-friendly correspondance between SAM flags and kinds of
# reads.
# The flags used below mean unmapped reads are excluded (-F 4) and only properly
# aligned pairs are kept (-f 3).
samtoolsReadFlags='-f 3 -F 4'

# See http://www.htslib.org/doc/samtools.html for a description of samtools
# mpileup options. Those used below mean that the minimum of the base quality
# the 'BAQ' quantity (see http://samtools.sourceforge.net/mpileup.shtml for an
# explanation) must be at least 5, and only the first 1000000 reads mapped to
# each point will be considered (NB a limit must be provided; the default is
# 250).
mpileupOptions='--min-BQ 5 --max-depth 1000000'

# Parameters for calling the consensus base at each position:
# The minimum coverage (number of reads) to call a base instead of a '?'
MinCov1=15
# The minimum coverage to use upper case for the base (to signal increased
# confidence)
MinCov2=30
# The minimum fraction of reads at a position before we call that base (or those
# bases, when one base alone does not reach that threshold fraction; e.g. say
# you have 60% A, 30% C and 10% G: if you set this fraction to 0.6 or lower we
# call an A, if you set it to 0.6-0.9 we call an M for "A or C", if you set it
# to 0.9-1 we call a V for "A, C or G".). Alternatively, if you choose a
# negative value, we always call the single most common base regardless of its
# fraction, unless two or more bases are equally (most) common, then we call the
# ambiguity code for those bases.
MinBaseFrac=-1

# Shall we remap to the consensus? (For remapping, gaps in coverage in the 
# consensus will filled in by the corresponding part of the orginal reference,
# and ambiguity codes will simplified to just one of the bases they represent.
# Because of this, if remapping to the consensus, you are strongly advised to 
# set the MinBaseFrac parameter above to any negative value.)
remap=true

# Shall we map contaminant reads to the reference (separately), to see which 
# reads would have contaminanted our final bam file had they not been removed?
MapContaminantReads=false

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
BestContigToRefAlignmentSuffix='_ContigsAndBestRef.fasta' # only for fully auto.
################################################################################
# The names of temporary files we'll create in the working directory.
# If you change the extension, you may well break something.
RawContigFile1='temp_HIVcontigs_uncut1.fasta'
RawContigFile2='temp_HIVcontigs_uncut2.fasta'
CutContigFile='temp_HIVcontigs_cut.fasta'
TempContigAlignment1='temp_HIVcontigs_wRefs_swap_1.fasta'
TempContigAlignment2='temp_HIVcontigs_wRefs_swap_2.fasta'
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
MapOutConversion1='temp_MapOutStep1'
MapOutConversion2='temp_MapOutStep2'
MapOutConversion3='temp_MapOutStep3'
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
ContigsWith1ref='temp_ContigsWith1ref.fasta'
RefMatchLog='temp_RefMatchLog.txt'
ContigAlignmentsToRefsDir='temp_ContigAlignmentsToRefsDir' # no whitespace!
SamtoolsSortFile='temp_SamtoolsSortFile'
