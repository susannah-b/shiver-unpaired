# shiver
Sequences from HIV Easily Reconstructed.
(Works for other viruses too.)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/shiver/README.html)  

<p align="center"><img src="docs/AssemblyPipelineDiagram_ForPaper.png" width=800, height="370"/></p>

shiver is a tool for assembling viral genomes by mapping reads (fragments of genetic sequence) to a custom reference sequence constructed using contigs specifc to each sample.
This reduces the biased loss of information that occurs from mapping to a reference that differs from the sample.
From the mapped reads, base frequencies are calculated, and a consensus sequence is called.
The method, its performance and scientific context are discussed [here](https://doi.org/10.1093/ve/vey007); please cite this if you find shiver, or our explanation of the relevant issues, helpful.
If you use shiver, please cite the paper for shiver and the papers of shiver's dependencies: details [here](docs/CitationDetails.bib).

shiver runs natively on Linux and Mac OS, but not Windows.

The easiest way of installing/updating shiver and all its dependencies is through [bioconda](https://bioconda.github.io/recipes/shiver/README.html).
Otherwise, after downloading the code from this repository (e.g. with `git clone https://github.com/ChrisHIV/shiver.git`), manually install all shiver's dependencies: [samtools](http://www.htslib.org/), [biopython](http://biopython.org/wiki/Download), [mafft](http://mafft.cbrc.jp/alignment/software/), [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) version 2.2.28 or higher (warning: earlier versions of blast have a bug that prevents shiver from correcting contigs), [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (optional: needed if you want to trim reads for quality or adapter sequences), [Fastaq](https://github.com/sanger-pathogens/Fastaq) (optional: needed if you want to trim reads for primers), and at least one of [smalt](http://www.sanger.ac.uk/science/tools/smalt-0) or [BWA](http://bio-bwa.sourceforge.net/) or [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) for mapping.
Installation instructions for all of these are [here](docs/InstallationNotes.sh).  

On any operating system (including Windows), if you have [VirtualBox](https://www.virtualbox.org/wiki/Downloads) installed, you can run [this](https://drive.google.com/file/d/1MohIOgJxcFVRv9v0aG2vyvEVgtm8lcyc/view?usp=sharing) image of Ubuntu Linux 22.04 which contains shiver and our separate tool [phyloscanner](https://github.com/BDI-pathogens/phyloscanner) (which allows you to investigate the within- and between-host diversity in mapped reads produced e.g. by shiver).
Instructions for the image (and full steps for how everything was installed) are [here](https://docs.google.com/document/d/1tX1juaiBhEZ5aW740ME05zaMm480Aj1eFvyY-H3Aor8/edit?usp=sharing).  

The shiver manual is [here](docs/ShiverManual.pdf).

shiver is freely available under the GNU General Public License version 3, described [here](LICENSE).  
