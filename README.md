# shiver
Sequences from HIV Easily Reconstructed.  

<p align="center"><img src="info/AssemblyPipelineDiagram_ForPaper.png" width=800, height="370"/></p>

shiver is freely available under the GNU General Public License version 3, described [here](LICENSE).  
shiver is a tool for mapping paired-end short reads to a custom reference sequence constructed using do novo assembled contigs, in order to minimise the biased loss of information that occurs from mapping to an reference that differs from the sample.
From the mapped reads, base frequencies are quantified, and a consensus sequence is called.  
The method and its performance are discussed [here](http://biorxiv.org/content/early/2016/12/09/092916); please cite this if you find shiver helpful.
If you use shiver, please also cite the publications of its dependencies. Citation details [here](info/CitationDetails.bib).

Dependencies: [Fastaq](https://github.com/sanger-pathogens/Fastaq), [samtools](http://www.htslib.org/), [biopython](http://biopython.org/wiki/Download), [mafft](http://mafft.cbrc.jp/alignment/software/), [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) version 2.2.28 or higher (warning: earlier versions of blast have a bug that prevents shiver from correcting contigs), [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (optional: needed if you want to trim reads for quality or adapter sequences) and at least one of [smalt](http://www.sanger.ac.uk/science/tools/smalt-0) or [BWA](http://bio-bwa.sourceforge.net/) or [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) for mapping.
Installation instructions for all of these are [here](info/InstallationNotes.sh).  

The shiver manual is [here](info/ShiverManual.pdf).
