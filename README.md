# shiver
Sequences from HIV Easily Reconstructed.  

<p align="center"><img src="info/AssemblyPipelineDiagram_ForPaper.png" width=800" height="370"/></p>

Dependencies: [smalt](http://www.sanger.ac.uk/science/tools/smalt-0), [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [Fastaq](https://github.com/sanger-pathogens/Fastaq), [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [samtools](http://www.htslib.org/), [biopython](http://biopython.org/wiki/Download), [mafft](http://mafft.cbrc.jp/alignment/software/).  
Before you begin processing a collection of samples there's an initialisation step: it should be run once only (i.e. not once for each sample).
It requires  
1. your choice of pipeline parameters, specified in `config.bash` (see section *The config file* later);  
2. an alignment of existing reference genomes (lots of which are available to download from the [Los Alamos National Lab](http://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html)) called `RefAlignment.fasta`, say;  
3. fasta files containing the adapters and primers used for sequencing (ask your sequencing team for these), called `adapters.fasta` and `primers.fasta`, say.
(Adapters are removed with trimmomatic, and "The naming of the various sequences within this file determines how they are used" - trimmomatic docs.
A default Illumina adapters file can be found in the [source code for IVA](https://github.com/sanger-pathogens/iva/).)  
Initialisation files will be put into a directory called `MyInitDir` if you run
```bash
$ ./shiver_init.bash MyInitDir config.bash RefAlignment.fasta adapters.fasta primers.fasta
```
(ignoring that first $ character, which just indicates that what follows should be run from command line). You shouldn't touch `MyInitDir` or the files therein after running this command, or unknown bad things might happen.

After running that initialisation command (once only!) you can process any number of samples. Since sample processing generates (many) files in the working directory, to avoid overwriting files you should work in an empty directory for each sample, as illustrated in the *Processing example 2* below.

#### Processing example 1: one sample, for which the contigs don't need correcting
Say you have forward reads in `reads_1.fastq.gz`, reverse reads in `reads_2.fastq.gz`, and contigs in `contigs.fasta`. Processing is achieved with the following two commands, replacing `SID` (sample ID) by a name used for labelling all output files for this sample:
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
and you're done. You'll probably want to delete the intermediate files produced during processing (beginning with 'temp' by default), i.e. `rm temp*`.

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

#### The config file
In `config.bash` you can change pipeline parameters from their default values. During development of shiver, `config.bash` has sometimes been updated (e.g. new features might require new parameters). In case this happens again in the future, to make sure you keep your changes it would be wise to make a copy of this file, change the copy, and use that instead of `config.bash` for your shiver commands. Then if `config.bash` changes in a future version of shiver, you can compare differences between the new version and your tweaked old version (e.g. with the `diff` command), and apply your changes to the new version.

#### Processing example 2: scripting over samples, all with HIV contigs, some of which need correcting
Before we start, note that scripted use of shiver (like any other command-line program) to process multiple files is easier if you know the basics of playing with filenames from the command line. Consider this toy example:
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

Now back to shiver. Say your directory looks like this:
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
$ for ContigFile in contigs/*.fasta; do  
  # Extract the SID from the filename: remove the extension and the path
  SID=$(basename "${ConfigFile%.fasta}")
  # Make, and change into, an empty directory for processing this contig file
  mkdir AlignmentOutput_"$SID"  
  cd AlignmentOutput_"$SID"  
  ../ShiverCode/shiver_align_contigs.bash ../MyInitDir ../config.bash "$ContigFile" "$SID"
  cd ..
done
```
For samples for which contig correction is necessary, `SID_cut_wRefs.fasta` will be produced as well as `SID_raw_wRefs.fasta` (the former containing the modified contigs).
In such cases, spliced contigs (those concatenating separated parts of the genome, detected by multiple separated blast hits) have been cut into pieces, and whole contigs or bits of contigs that are reverse complements (detected by inverted blast hits) are reverse complemented back to normal. Our detection of contigs needing correction does not have 100% specificity - e.g. a large indel in a contig can cause two separate blast hits - and so where `SID_cut_wRefs.fasta` exists you're advised to inspect both it and `SID_raw_wRefs.fasta` (instead of blindly trusting that the cut one is better), choose the better looking one and discard the other, and then delete ragged ends if they appear (as previously). Let's say you put all these files together in a directory called `CheckedContigAlignments`.  
```bash
# For samples that had an SID_cut_wRefs.fasta file, we kept either that file or
# SID_raw_wRefs.fasta. Let's rename all these files to have the same format,
# removing '_raw' or '_cut' to leave just SID_wRefs.fasta, to make it easier
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
  ../ShiverCode/shiver_map_reads.bash ../MyInitDir ../config.bash "$ContigFile" "$SID" "$BlastFile" "$alignment" "$reads1" "$reads2"
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
For example in processing example 2, after `alignment=../CheckedContigAlignments/"$SID"_wRefs.fasta`, you could have
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
Decreasing these thresholds increases sensitivity to HIV reads at the cost of specificity (since some reads are bound to be contaminants), and you ought to balance these, not just maximimise the former.
For example with dated sequences, you could try to get the strongest correlation between real time and the evolutionary distance inferred from a phylogeny - the R^2 of the molecular clock - since too little real sequence and too much contaminant sequence will both screw up your phylogeny.


##### TODO: say something about global alignment coords
