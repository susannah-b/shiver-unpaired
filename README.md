# shiver
Sequences from HIV Easily Reconstructed.  

<p align="center"><img src="info/AssemblyPipelineDiagram_ForPaper.png" width=800, height="370"/></p>

Dependencies: [smalt](http://www.sanger.ac.uk/science/tools/smalt-0), [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [Fastaq](https://github.com/sanger-pathogens/Fastaq), [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [samtools](http://www.htslib.org/), [biopython](http://biopython.org/wiki/Download), [mafft](http://mafft.cbrc.jp/alignment/software/).  
Before you begin processing a collection of samples there's an initialisation step: it should be run once only (i.e. not once for each sample).
It requires  
1. your choice of pipeline parameters, specified in `config.sh` (see section *The config file* later);  
2. an alignment of existing reference genomes (lots of which are available to download from the [Los Alamos National Lab](http://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html)) called `RefAlignment.fasta`, say;  
3. fasta files containing the adapters and primers used for sequencing (ask your sequencing team for these), called `adapters.fasta` and `primers.fasta`, say.
(Adapters are removed with trimmomatic, and "The naming of the various sequences within this file determines how they are used" - trimmomatic docs.
A default Illumina adapters file can be found in the [source code for IVA](https://github.com/sanger-pathogens/iva/).)  
Initialisation files will be put into a directory called `MyInitDir` if you run
```bash
$ ./shiver_init.sh MyInitDir config.sh RefAlignment.fasta adapters.fasta primers.fasta
```
(ignoring that first $ character, which just indicates that what follows should be run from command line).
You shouldn't touch `MyInitDir` or the files therein after running this command, or unknown bad things might happen.

After running that initialisation command (once only!) you can process any number of samples.
Since sample processing generates (many) files in the working directory, to avoid overwriting files you should work in an empty directory for each sample, as illustrated in the *Scripted usage to process batches of samples* section below.

#### How to process a single sample
Say you have forward reads in `reads_1.fastq.gz`, reverse reads in `reads_2.fastq.gz`, and contigs in `contigs.fasta`.
Processing requires two commands.
The first is  
```bash
$ ./shiver_align_contigs.sh MyInitDir config.sh contigs.fasta SID
```
replacing `SID` (sample ID) by a name used for labelling all output files for this sample.
This command will produce a file named `SID.blast` (detailing blast hits of your contigs to those existing references supplied for the initialisation).
Assuming at least one contig looks like HIV, `SID.blast` will not be empty (see section *Including samples without contigs* otherwise), and there will be another file called `SID_raw_wRefs.fasta` - an alignment of the HIV contigs (i.e. those that blasted) to your input existing reference genomes.
As discussed in the shiver article, this alignment should be inspected, and in the minority of cases where this is required, edited; see the following section.

If a single contig has two or more blast hits (ignoring hits wholly inside other hits), it is taken to be chimeric / spliced: erroneously connecting disconnected parts of the genome.
shiver cuts such contigs between the hits, generating two or more contigs that can be aligned separately, and so for such samples there is an additional file `SID_cut_wRefs.fasta` containing the cut contigs (and any that did not need to be cut) aligned to the existing references.
`SID_cut_wRefs.fasta` *should not be assumed to be superior to* `SID_raw_wRefs.fasta`.
For example unusual indels can result in multiple blast hits, as can close matches to shorter references combined with less close matches to longer references.
You should therefore look at both alignments, choose which set of contigs you think are a closer approximation to the sample consensus, edit if needed, and discard the other. 

With an alignment of contigs you're happy with, say `SID_raw_wRefs.fasta`, run
```bash
$ ./shiver_map_reads.sh MyInitDir config.sh contigs.fasta SID SID.blast SID_raw_wRefs.fasta reads_1.fastq.gz reads_2.fastq.gz
```
and you're done.
You'll probably want to delete the intermediate files produced during processing (beginning with 'temp' by default), i.e. `rm temp*`.

#### A manual step between the two automatic steps? Blasphemy!
Correct alignment has already been attempted automatically in the first step, but unfortunately aligning HIV whole genomes automatically, perfectly, 100% of the time is a dream: HIV has one of the highest mutation rates known to biology, including frequent insertions and deletions.
Compounding this problem is the fact that de novo assembly output is sometimes an imperfect representation of the sample, and aligning problematic sequences may give problematic alignments. 
For most samples the contigs and their alignment are fine, and seeing this from the alignment takes one or two seconds of human time.
When editing is required, it amounts to one of two things: deleting a whole contig if it appears to be wholly junk (i.e. short stretches of poorly aligned sequence separated by large gaps), or deleting the end of a contig if just the end looks wrong (i.e. a stretch of sequence in very poor agreement with the existing references, or a stretch of sequence separated from the rest of the contig by a gap implausibly large for a real deletion).
Visually checking and editing when needed ensures an alignment of contigs you can trust, allowing the second `shiver` command to reliably construct a) a tailored reference that minimises mapping bias, and b) a global alignment of all samples, instantly, without further need for an alignment algorithm: see the end of section *Scripted usage to process batches of samples*.  
[This](https://github.com/olli0601/PANGEAhaircut) R package by Oliver Ratmann uses machine learning to correct such alignments; however for each sample you must manually check the corrections to ensure they are correct, instead of checking the alignment itself.

#### What output do I get?

* `SID.bam`, `SID_ref.fasta`: the bam file of mapped reads, and the reference to which they're mapped (it's a shame the bam file format does not include the reference but that's how it is).
* `SID_BaseFreqs.csv`: the frequencies of A, C, G, T, '-' for a deletion, and N for unknown, at each position in the genome.
* `SID_MinCov_foo_bar.fasta`: the pairwise alignment of the consensus genome and the reference used for mapping.
In place of `foo` and `bar` here you'll see the two coverage thresholds specified in `config.sh`: `foo` is the minimum number of reads to call a base (if there are fewer than `foo` reads we call `?` instead of a base), `bar` the minimum number to use upper case for the base (to signal increased confidence).
The mapping reference is included with the consensus here to give context to any `?` characters.
* `SID_MinCov_foo_bar_wContigs.fasta`: an alignment of the consensus, the reference used for mapping, and all HIV contigs for this sample (as they were after any automatic cutting and any manual editing).
* `SID_MinCov_foo_bar_ForGlobalAln.fasta`, `SID_coords.csv`: these files are useful for when multiple samples are processed, so we postpone their explanation.  
* `SID_clean_1.fastq`, `SID_clean_2.fastq`: the unmapped reads after removal of adapters, primers, low-quality bases, and those read pairs suspected of being contamination.  
* `SID_InsertSizeCounts.csv`: the insert-size distribution.
* `SID_ContaminantReads.bam`: a bam file of only the contaminant reads mapped to the HIV reference (to illustrate the importance of removing these reads from the HIV reads prior to mapping).

#### The config file
In `config.sh` you can change pipeline parameters from their default values.
During development of shiver, `config.sh` has sometimes been updated (e.g. new features might require new parameters).
In case this happens again in the future, to make sure you keep your changes it would be wise to make a copy of this file, change the copy, and use that instead of `config.sh` for your shiver commands.
Then if `config.sh` changes in a future version of shiver, you can compare differences between the new version and your tweaked old version (e.g. with the `diff` command), and apply your changes to the new version.

#### Scripted usage to process batches of samples
Before we start, note that scripted use of shiver (like any other command-line program) to process multiple files is easier if you know the basics of playing with filenames from the command line.
Consider this toy example:
```bash
# We have two directories: one with text files, one with csv files.
$ ls CsvDir/
bar.csv  foo.csv
$ ls TxtDir/
bar.txt  foo.txt
# Assign the text file we want to a variable:
$ DesiredTextFile=TxtDir/foo.txt
# This shows how to find the associated csv file:
$ DesiredTextFileNoPath=$(basename "$DesiredTextFile")
$ echo "$DesiredTextFileNoPath"
foo.txt
$ DesiredTextFileNoPathNoExtention="${DesiredTextFileNoPath%.txt}"
$ echo "$DesiredTextFileNoPathNoExtention"
foo
$ AssociatedCsvFile=CsvDir/"$DesiredTextFileNoPathNoExtention".csv
$ ls "$AssociatedCsvFile"
CsvDir/foo.csv
# Or, the above steps all in one go:
$ DesiredTextFile=TxtDir/foo.txt
$ ls CsvDir/$(basename "${DesiredTextFile%.txt}").csv
CsvDir/foo.csv
```

Now back to shiver.
Say your directory looks like this:
```bash
$ ls
contigs
MyInitDir
reads
ShiverCode
config.sh
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
$ for ContigFile in contigs/*.fasta; do  
  # Extract the SID from the filename: remove the extension and the path
  SID=$(basename "${ConfigFile%.fasta}")
  # Make, and change into, an empty directory for processing this contig file
  mkdir AlignmentOutput_"$SID"  
  cd AlignmentOutput_"$SID"  
  ../ShiverCode/shiver_align_contigs.sh ../MyInitDir ../config.sh "$ContigFile" "$SID"
  cd ..
done
```
As explained above, for those samples for which contig correction is necessary, `SID_cut_wRefs.fasta` will be produced as well as `SID_raw_wRefs.fasta`, and only the better looking of these two should be kept (and edited if needed).
Let's say you put one of these two files for each sample all together in a directory called `CheckedContigAlignments`.  
```bash
# For samples that had an SID_cut_wRefs.fasta file, we kept either that file or
# SID_raw_wRefs.fasta. Let's rename all these files to have the same suffix -
# removing '_raw' or '_cut' to leave just SID_wRefs.fasta - to make it easier
# to find the one we want.  
$ cd CheckedContigAlignments  
$ for alignment in *_cut_wRefs.fasta; do
  mv -i "$alignment" "${alignment%_cut_wRefs.fasta}"_wRefs.fasta
done
$ for alignment in *_raw_wRefs.fasta; do
  mv -i "$alignment" "${alignment%_raw_wRefs.fasta}"_wRefs.fasta
done
# Now let's map! Make and change into an empty directory for processing each sample.
$ cd ..
$ for ContigFile in contigs/*.fasta; do  
  # Find the other files for this sample.
  # NB careful scripting would check that files exist before using them.
  SID=$(basename "${ConfigFile%.fasta}")
  mkdir MappingOutput_"$SID"
  cd MappingOutput_"$SID"
  BlastFile=../AlignmentOutput_"$SID"/"$SID".blast
  alignment=../CheckedContigAlignments/"$SID"_wRefs.fasta
  reads1=../reads/"$SID"_1.fastq.gz
  reads2=../reads/"$SID"_2.fastq.gz
  ../ShiverCode/shiver_map_reads.sh ../MyInitDir ../config.sh "$ContigFile" "$SID" "$BlastFile" "$alignment" "$reads1" "$reads2"
  cd
done
```
and you're done mapping.

Now for those `*_MinCov_foo_bar_ForGlobalAln.fasta` files.
These are generated by excising unique insertions seen in a sample (i.e. not seen in any of the existing references used for initialisation), and then inserting gaps using coordinate translation in such a way that these files are all aligned with each other and with the alignment of existing references.
Constructing a global alignment is then achieved just by putting all those files into one file:
```bash
$ cat MappingOutput_*/*_ForGlobalAln.fasta > GlobalAln.fasta
$ cat MyInitDir/ExistingRefAlignment.fasta >> GlobalAln.fasta
```
Now it's over to you for phylogenetics, GWAS etc.

#### Including samples without contigs
You might have some samples for which none of the contigs blast to the existing references you supplied, i.e. it looks like all contigs are contaminants.
In this case the contig alignment step will produce an empty blast file and no `SID_raw_wRefs.fasta` file.
You could choose to ignore such cases - that's reasonable.
You could also choose to try to recover some HIV reads in your contamination-dominated sample.
To this end, in place of the alignment of contigs to existing references amongst the arguments to `shiver_map_reads.sh`, you can supply a fasta file containing a single sequence: that sequence will be used as the reference for mapping, instead of a tailored one constructed from contigs.
If that sequence is one of the existing references you provided at the initialisation step, `shiver` knows how to do the coordinate translation necessary to produce a `SID_MinCov_foo_bar_ForGlobalAln.fasta` file; if not, this file will not be produced (you'll still get your consensus sequence though).

To script this kind of thing, you can just check whether the alignment of contigs to existing references exists for this sample: if not, choose a sequence to use as your reference to mapping.
This would probably be most easily achieved by making a big look-up table before you start.
(e.g. running the program [kraken](https://ccb.jhu.edu/software/kraken/) on the reads for each sample, one can see which of the existing references has the largest number of reads attributed to it; one might obtain a non-null result even when no HIV contigs were assembled.)
For example in the scripted usage above, after `alignment=../CheckedContigAlignments/"$SID"_wRefs.fasta`, you could have
```bash
  if [ ! -f "$alignment" ]; then
    # Reassign alignment to be the reference previously chosen for this sample
  fi
```

#### Analysis and partial reprocessing

##### Using different coverage thresholds

Recall `foo` and `bar` above.
You can generate another consensus sequence using different values of `foo` and `bar` without rerunning the whole pipeline thus:
```bash
$ tools/CallConsensus.py SID_BaseFreqs.csv foo bar > MyNewConsensus.fasta
```
If you decrease the value of `foo` you will see more bases and fewer `?` characters: the less strict you are with your requirement on 'vertical' coverage (the height at a given position in the bam file, if you like), the more 'horizontal' coverage (i.e. genomic coverage, the total number of bases called) will increase.
That's good right?
Yes and no.
You could also increase your horizontal coverage by replacing all `?` characters by a random selection of As, Cs, Gs and Ts, but it's not recommended.
Decreasing these thresholds increases sensitivity to HIV reads at the cost of specificity (since some reads are bound to be contaminants), and you ought to balance these, not just maximimise the former.
For example with dated sequences, you could try to get the strongest correlation between real time and the evolutionary distance inferred from a phylogeny - the R^2 of the molecular clock - since too little real sequence and too much contaminant sequence will both screw up your phylogeny.
To regenerate a coordinate-translated version of this consensus for the global alignnent (of all consensuses produced by shiver), `tools/TranslateSeqForGlobalAln.py` can be run, taking as its two arguments the consensus, and `SID_coords.csv` generated by the full run of shiver.

##### Discarding dissimilar reads

One parameter in config.sh is the minimum read identity – the fraction of bases in the read which
are mapped and agree with the reference – required for a read to be considered mapped, and so retained in the bam file.
If you wish to increase this after completion of shiver, reads with an identity below your new higher threshold can be discarded by running `tools/RemoveDivergentReads.py` on a bam file.
Running `shiver_reprocess_bam.sh` on the resulting bam file (or indeed any bam file) implements just the last steps in shiver, namely generating pileup, calculating the base frequencies, and calling the consensus.

##### Number and accuracy of mapped reads

`tools/FindNumMappedBases.py` calculates the total number of mapped nucleotides in a bam file (i.e. the number of mapped reads multiplied by read length, minus the total length of sequence clipped from reads) optionally binned by read identity.
In the absence of mapped contaminant reads, and all else being equal, mapping to a reference which is closer to the true consensus should map more nucleotides and mapped reads should have higher identities.

##### Finding mapping problems

`tools/FindClippingHotSpots.py` counts, at each position in the genome, the number and percentage of reads that are clipped from that position to their left or right end.
(Optionally with a minimum length of clipped sequence for it to be counted.)
Having many such reads is a warning sign of a structural difference (an indel) between the reference and the reads which is hindering the latter from being mapped accurately to the former.
