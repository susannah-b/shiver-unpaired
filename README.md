# shiver
Sequences from HIV Easily Reconstructed.  

<p align="center"><img src="info/AssemblyPipelineDiagram_ForPaper.png" width=800" height="370"/></p>

Dependencies: [smalt](http://www.sanger.ac.uk/science/tools/smalt-0), [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [Fastaq](https://github.com/sanger-pathogens/Fastaq), [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) [samtools](http://www.htslib.org/), [biopython](http://biopython.org/wiki/Download), [mafft](http://mafft.cbrc.jp/alignment/software/).  
Before you begin processing samples there's a one-off initialisation step.
This requires  
1. your choice of pipeline parameters specified in `config.bash`,  
2. an alignment of existing reference genomes, lots of which are available to download from the [Los Alamos National Lab](http://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html),  
3. fasta files containing the adapters and primers used for sequencing: ask your sequencing team for these.  
Initialisation files will be put into a directory called `MyInitDir` if you run
```bash
$ ./shiver_init.bash MyInitDir config.bash MyRefAlignment.fasta
```
Sample processing generates (many) files in the working directory, so to avoid overwriting files it's recommended to work in an empty directory.

#### Processing example 1: one sample, for which the contigs don't need correcting
Say you have forward reads in `reads_1.fastq.gz`, reverse reads in `reads_1.fastq.gz`, and contigs in `contigs.fasta`. Processing is achieved with the following two commands, replacing `SID` (sample ID) by name for labelling all your output files:
```bash
$ ./shiver_align_contigs.bash MyInitDir config.bash contigs.fasta SID
```
Amongst the files this step produces is `SID_raw_wRefs.fasta` - an alignment of your HIV contigs to your input existing reference genomes - which you should should visually check, deleting ragged ends which arise occasionally.
Also produced is `SID.blast` (the result of blasting your contigs to those existing references.
Now run
```bash
$ ./shiver_map_reads.bash MyInitDir config.bash contigs.fasta SID SID.blast SID_raw_wRefs.fasta reads_1.fastq.gz reads_2.fastq.gz
```
and you're done.

#### What output do I get?

`SID.bam`, `SID_ref.fasta`: the bam file of mapped reads, and the (tailored) reference to which they're mapped.  
`SID_BaseFreqs.csv`: the frequencies of the four bases (plus the `-` character for deletions) at each position in the bam file.  
`SID_MinCov_foo_bar.fasta`: the consensus genome. `foo` and `bar` here are the two coverage thresholds specified in `config.bash`: `foo` is the minimum number of reads to call a base, `bar` the minimum number to use upper case for the base (an indicator of increased confidence).  
`SID_MinCov_foo_bar_ForGlobalAln.fasta`, `SID_coords.csv`: these files are useful for when multiple samples are processed, so we postpone their explanation.
`SID_clean_1.fastq.gz`, `SID_clean_2.fastq.gz`: the reads after removal of adapters, primers, low-quality bases, and those read pairs suspected of being contamination.

#### Processing example 2: scripting over samples, all with HIV contigs, some of which need correcting
Say your directory looks like this:
```bash
$ ls
contigs
MyInitDir
reads
ShiverCode
config.bash
```
with the first four being subdirectories containing what you'd expect.
Here are your contigs:
```bash
$ ls contigs/
sampleA.fasta sampleB.fasta sampleC.fasta
```
and here are your reads:
```bash
$ ls reads/
sampleA_1_fastq.gz sampleA_2_fastq.gz  
sampleB_1_fastq.gz sampleB_2_fastq.gz  
sampleC_1_fastq.gz sampleC_2_fastq.gz
```
Start by aligning the contigs for all samples:
```bash
$ mkdir AlignmentOutput  
$ cd AlignmentOutput  
$ for ContigFile in ../contigs/*.fasta; do  
  # Extract the SID from the filename: remove the extension and the path
  SID=$(basename "${ConfigFile%.fasta}")
  ../ShiverCode/shiver_align_contigs.bash ../MyInitDir ../config.bash "$ContigFile" "$SID"
done
```
For samples for which contig correction is necessary, `SID_cut_wRefs.fasta` will be produced as well as `SID_raw_wRefs.fasta` (the former containing the modified contigs).  
In such cases, spliced contigs (those concatenating separed parts of the genome, detected by multiple separated blast hits) have been cut into pieces, and reverse-complement contigs (detected by inverted blast hits) are reverse complemented.
This doesn't work perfectly every time -- a large indel in a contig can cause two separate blast hits, and cutting the contig at slightly the wrong spot can cause alignment difficulty -- and so where `SID_cut_wRefs.fasta` exists you're advised to inspect both it and `SID_raw_wRefs.fasta`.

Running the aligning step,

Files that definitely will be produced:
SID.blast

Files that might be produced:
SID_hiv.fasta, SID_raw_wRefs.fasta: if there's at least one HIV contig
SID_hiv_cut.fasta, SID_cut_wRefs.fasta: if contig correction was necessary

Running the mapping step,

Files that definitely will be produced:
SID.bam, SID.bam.bai
SID_ref.fasta, SID_ref.fasta.fai
SID_BaseFreqs.csv
SID_InsertSizeCounts.csv

SID_MinCov_foo_bar.fasta, where foo and bar are your chosen thresholds on coverage to call a base at all, and to call a base with upper case, respectively.
SID_ContaminantReads.bam: if there were no contaminant contigs, or if no read pair blasted to the contaminant contigs better than to the tailored reference, the file will still be produced but as a blank file. This is to make it more easy to find the distribution over all your samples of the number of contaminant reads: you just count the number of reads present in each of these files.

Files that might be produced:
SID_clean_1.fastq.gz, SID_clean_2.fastq.gz: if there was at least one contaminant contig and at least one read pair blasted better to it/them than to the tailored reference.
SID_MinCov_foo_bar_ForGlobanAln.fasta, SID_coords.csv: only won't be produced if you manually specify a reference for mapping (instead of constructing one, e.g. if there are no HIV contigs) and you chose one which is not in your input alignment of existing references. Again, foo and bar are your chosen thresholds on coverage.
