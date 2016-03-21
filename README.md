# shiver
Sequences from HIV Easily Reconstructed.  

<p align="center"><img src="info/AssemblyPipelineDiagram_ForPaper.png" width=800" height="370"/></p>

Dependencies: [smalt](http://www.sanger.ac.uk/science/tools/smalt-0), [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [Fastaq](https://github.com/sanger-pathogens/Fastaq), [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [samtools](http://www.htslib.org/), [biopython](http://biopython.org/wiki/Download), [mafft](http://mafft.cbrc.jp/alignment/software/).  
Before you begin processing samples there's a one-off initialisation step.
This requires  
1. your choice of pipeline parameters, specified in `config.bash`;  
2. an alignment of existing reference genomes, lots of which are available to download from the [Los Alamos National Lab](http://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html);  
3. fasta files containing the adapters and primers used for sequencing: ask your sequencing team for these.
(Adapters are removed with trimmomatic, and "The naming of the various sequences within this file determines how they are used" - trimmomatic docs.
A default Illumina adapters file can be found in the [source code for IVA](https://github.com/sanger-pathogens/iva/).)  
Initialisation files will be put into a directory called `MyInitDir` if you run
```bash
$ ./shiver_init.bash MyInitDir config.bash MyRefAlignment.fasta
```
(ignoring that first $ character, which just indicates this should be run from command line).   Sample processing generates (many) files in the working directory, so to avoid overwriting files you're advised to work in an empty directory.

#### Processing example 1: one sample, for which the contigs don't need correcting
Say you have forward reads in `reads_1.fastq.gz`, reverse reads in `reads_2.fastq.gz`, and contigs in `contigs.fasta`. Processing is achieved with the following two commands, replacing `SID` (sample ID) by a name used for labelling all your output files:
```bash
$ ./shiver_align_contigs.bash MyInitDir config.bash contigs.fasta SID
```
Amongst the files this step produces is `SID_raw_wRefs.fasta` - an alignment of your HIV contigs to your input existing reference genomes - which you should should visually check, and manually delete ragged contig ends which arise occasionally.
(NB you may modify the contigs but not the existing reference sequences in the alignment.)
Also produced is `SID.blast` - the result of blasting your contigs to those existing references.
Now run
```bash
$ ./shiver_map_reads.bash MyInitDir config.bash contigs.fasta SID SID.blast SID_raw_wRefs.fasta reads_1.fastq.gz reads_2.fastq.gz
```
and you're done.

#### I notice there's a manual step between the two automatic steps. Blasphemy!
Correct alignment has already been attempted automatically in the first step, but unfortunately aligning HIV sequences automatically, perfectly, 100% of the time is a dream: HIV has the highest mutation rate known to biology, and the alignment algorithm has an unimaginably large space of possibilities to explore.
Most of the time it works, and inspecting the alignment it will take you one or two seconds to see that that's the case.
Most of the cases in which it fails, all you have to do is to delete a short stretch of sequence separated from the main body of the contig by a long erroneous gap, and possibly delete a whole contig if it scored a blast hit (and so was considered to be HIV, and was aligned) but is actually not HIV, so when 'aligned' is just a mess.
Performing this step manually ensures an alignment you can trust, allowing the second `shiver` command to reliably construct a) a tailored reference that minimises mapping bias, and b) a global alignment of all samples, instantly, without further need for an alignment algorithm: see the end of section *Processing example 2*.  
[This](https://github.com/olli0601/PANGEAhaircut) R package by Oliver Ratmann uses machine learning to correct these cases of failed alignment; however for each sample you must manually check the corrections to ensure they are correct, instead of checking the alignment itself.



#### What output do I get?

`SID.bam`, `SID_ref.fasta`: the bam file of mapped reads, and the (tailored) reference to which they're mapped.  
`SID_BaseFreqs.csv`: the frequencies of the four bases (plus the `-` character for deletions) at each position in the bam file.  
`SID_MinCov_foo_bar.fasta`: the consensus genome.
In place of `foo` and `bar` here you'll see the two coverage thresholds specified in `config.bash`: `foo` is the minimum number of reads to call a base (if there are fewer than `foo` reads we call `?` instead of a base), `bar` the minimum number to use upper case for the base (to signal increased confidence).  
`SID_MinCov_foo_bar_ForGlobalAln.fasta`, `SID_coords.csv`: these files are useful for when multiple samples are processed, so we postpone their explanation.  
`SID_clean_1.fastq.gz`, `SID_clean_2.fastq.gz`: the reads after removal of adapters, primers, low-quality bases, and those read pairs suspected of being contamination.  
`SID_InsertSizeCounts.csv`: the inferred insert-size distribution.

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
with the first four being subdirectories containing what the name suggests.
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
In such cases, spliced contigs (those concatenating separated parts of the genome, detected by multiple separated blast hits) have been cut into pieces, and reverse-complement contigs (detected by inverted blast hits) are reverse complemented.
This doesn't work perfectly every time - a large indel in a contig can cause two separate blast hits, and cutting the contig at slightly the wrong spot can cause alignment difficulty - and so where `SID_cut_wRefs.fasta` exists you're advised to inspect both it and `SID_raw_wRefs.fasta`, choose the better looking and discard the other, and delete ragged ends if they appear (as previously).
```bash
# For samples that had an SID_cut_wRefs.fasta file, we kept either that file or
# SID_raw_wRefs.fasta. Let's rename all these files so they have the same format,
# namely SID_wRefs.fasta.
$ for alignment in *_cut_wRefs.fasta; do
  mv -i "$alignment" "${alignment%_cut_wRefs.fasta}"_wRefs.fasta
done
$ for alignment in *_raw_wRefs.fasta; do
  mv -i "$alignment" "${alignment%_raw_wRefs.fasta}"_wRefs.fasta
done
# Now let's map!
$ cd ..
$ mkdir MappingOutput
$ cd MappingOutput
$ for ContigFile in ../contigs/*.fasta; do  
  # Find the other files for this sample.
  # NB careful scripting would check that files exist before using them.
  SID=$(basename "${ConfigFile%.fasta}")
  BlastFile=../AlignmentOutput/"$SID".blast
  alignment=../AlignmentOutput/"$SID"_wRefs.fasta
  reads1=../reads/"$SID"_1.fastq.gz
  reads2=../reads/"$SID"_2.fastq.gz
  ../ShiverCode/shiver_map_reads.bash ../MyInitDir ../config.bash "$ContigFile" "$SID" "$BlastFile" "$alignment" "$reads1" "$reads2"
done
```
and you're done mapping.

Now for those `*_MinCov_foo_bar_ForGlobalAln.fasta` files.
These are generated by excising unique insertions seen in a sample (i.e. not seen in any of the existing references used for initialisation), and then inserting gaps using coordinated translation in such a way that these files are all aligned with each other and with the alignment of existing references.
Constructing a global alignment is then as easy as this:
```bash
$ cat *_ForGlobalAln.fasta > GlobalAln.fasta
$ cat ../MyInitDir/ExistingRefAlignment.fasta >> GlobalAln.fasta
```
Now it's over to you for phylogenetics etc.

#### Processing example 3: scripting when some samples have no HIV contigs
You might have some samples for which none of the contigs blast to the existing references you supplied, i.e. it looks like all contigs are contaminants.
In this case the contig alignment step will produce an empty blast file and no `SID_raw_wRefs.fasta` file.
You could choose to ignore such cases - that's reasonable.
You could also choose to try to recover some HIV reads in your contamination-dominated sample.
To this end, in place of the alignment of contigs to existing references amongst the arguments to `shiver_map_reads.bash`, you can supply a fasta file containing a single sequence: that sequence will be used as the reference for mapping, instead of a tailored one constructed from contigs.
If that sequence is one of the existing references you provided at the initialisation step, `shiver` knows how to do the coordinate translation necessary to produce a `SID_MinCov_foo_bar_ForGlobalAln.fasta` file; if not, this file will not be produced (you'll still get your consensus sequence though).

To script this kind of thing, you can just check whether the alignment of contigs to existing references exists for this sample: if not, choose a sequence to use as your reference to mapping.
This would probably be most easily achieved by making a big look-up table before you start.
(e.g. running the program [kraken](https://ccb.jhu.edu/software/kraken/) on the reads for each sample, one can see which of the existing references has the largest number of reads attributed to it; from experience one often obtains a non-null result even when no HIV contigs were assembled.)
For example in processing example 2, after `alignment=../AlignmentOutput/"$SID"_wRefs.fasta`, you could have
```bash
  if [ ! -f "$alignment" ]; then
    # Reassign alignment to be the reference previously chosen for this sample
  fi
```

#### Using different coverage thresholds
Recall `foo` and `bar` above.
You can generate another consensus sequence using different values of `foo` and `bar` without rerunning the whole pipeline thus:
```bash
$ tools/CallConsensus.py SID_BaseFreqs.csv foo bar > MyNewConsensus.fasta
```
If you decrease the value of `foo` you will see more bases and fewer `?` characters: the less strict you are with your requirement on 'vertical' coverage (the height at a given position in the bam file, if you like), the more 'horizontal' coverage (i.e. genomic coverage, the total number of bases called) will increase.
That's good right?
Yes and no.
You could also increase your horizontal coverage by replacing all `?` characters by a random selection of As, Cs, Gs and Ts, but it's not recommended.
Decreasing these thresholds increases sensitivity to HIV reads at the cost of specificity (since some reads are bound to be contaminants).
Choose your thresholds to maximise the information you get out of a _batch_ of samples, not a single sample.
(For example with dated sequences, you could try balancing sensitivity and specificity to get the strongest correlation between real time and the evolutionary distance inferred from a phylogeny - the R^2 of the molecular clock - since too little real sequence and too much contaminant sequence will both screw up your phylogeny.)


##### TODO: say something about global alignment coords
##### Some extra info, not properly formatted yet
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
SID_MinCov_foo_bar_ForGlobanAln.fasta, SID_coords.csv: only won't be produced if you manually specify a reference for mapping (instead of constructing one, e.g. if there are no HIV contigs) and you chose one which is not in your input alignment of existing references. Again, foo and bar are your chosen thresholds on coverage.
