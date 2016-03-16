#------------------------------------------------------------------------------------------------#
####                ANALYZE OUTPUT OF BLAST AND SAVE ONLY HIV SEQUENCES                       ####
#------------------------------------------------------------------------------------------------#

rm(list = ls())
# get arguments
args <- commandArgs(trailingOnly = TRUE)

blast.filename <- args[1] # "~/DID/BEEHIVE/Collective_notebook/blast.08.12.2014"
contigs.filename <- args[2] # "~/DID/BEEHIVE/alignments/all.contigs.08.12.2014"
alignment.out.name <- args[3] # "~/DID/BEEHIVE/alignments/test01.04.2015.fasta"

sort.blast.output <- function(blast.table, contigs){
  
  stopifnot(names(blast.table) == c("Contig_name","Match","Eval","PropI","QStart","QEnd","MStart","MEnd")) # check the table has the correct columns
  
  # this function uses the blast output and returns a list of cleaned contigs
  # cleaning includes elimination of short blast matches,
  # reverse - complementation of contigs when necesary
  # cutting and reverse-complementation of contigs when necesary
  
  # remove all matches <= 100 bp
  blast.table <- blast.table[which(blast.table$QEnd - blast.table$QStart > 100), ]
  # TODO: exit if there are no matches left after this step.
  
  all.contigs.id <- unique(blast.table$Contig_name) # name of all contigs
  
  new.contigs <- list() # list of sequences of contigs once they are sorted
  names.new.contigs <- c() # list of names
  n.contigs <- 0
  
  if(!all(all.contigs.id %in% names(contigs)))stop("some contigs are present in blast table but not present in contigs")
  
  for(ii in 1:length(all.contigs.id)){# for each contig
    
    idx <- which(blast.table$Contig_name == all.contigs.id[ii]) # idx of contigs in blast.matches
    idx.in.seq <- which(names(contigs) == all.contigs.id[ii]) # idx of contigs in the list of sequences
    sub.tab <- blast.table[idx,] # the corresponding blast outputs
    if(nrow(sub.tab)>1)sub.tab <- remove.subset(sub.tab) # eliminate the blast outputs which are a subset of another
    tot.size.region <- tot.match.length(sub.tab) # calculate total size of region covered by the matches
    sub.tab <- sub.tab[order(sub.tab$QStart),]
    
    if(nrow(sub.tab) == 1){ # one match only: easy - just add the contig if it is long enough, reverse if necesary
      
      if(tot.size.region > 0.7 * length(contigs[[idx.in.seq]])){ # if the match is big enough
        n.contigs = n.contigs+1
        if(sub.tab$MStart[1] > sub.tab$MEnd[1]){ # the match is in reverse direction
          new.contigs[[n.contigs]] <- reverse.complement(contigs[[idx.in.seq]])
          cat("contig ", as.character(all.contigs.id[ii]), " reversed\n")
        } else { # the match is in correct direction
          new.contigs[[n.contigs]] <- contigs[[idx.in.seq]]
        }        
        names.new.contigs <- c(names.new.contigs, names(contigs)[idx.in.seq])
      }
      
    } else { # several matches: cut!
      
      cat("contig ", as.character(all.contigs.id[ii]), " cut in", nrow(sub.tab), "\n")
      cutpoints <- c()
      for(kk in 1:(nrow(sub.tab)-1)){
        cutpoints <- c(cutpoints, floor(0.5*(sub.tab[kk,]$QEnd + sub.tab[kk+1,]$QStart)))
      }
      for(kk in 1:nrow(sub.tab)){ # add cut contigs one by one
        if(kk == 1) begin = 1 else begin = cutpoints[kk-1] + 1
        if(kk == nrow(sub.tab)) end = length(contigs[[idx.in.seq]]) else end = cutpoints[kk]
        n.contigs = n.contigs+1
        if(sub.tab$MStart[kk] < sub.tab$MEnd[kk]){ # this piece does not need inversion
          new.contigs[[n.contigs]] <- contigs[[idx.in.seq]][begin:end]
        } else {# this piece needs inversion
          new.contigs[[n.contigs]] <- reverse.complement(contigs[[idx.in.seq]][begin:end])
        }
        names.new.contigs <- c(names.new.contigs, paste(names(contigs)[idx.in.seq], ".", kk, sep = ""))
      }    
    }  
  }
  return(list(new.contigs, names.new.contigs))
}

tot.match.length <- function(sub.table){ # returns the total sequence length that is covered by blast
  set.pos <- sub.table$QStart[1]:sub.table$QEnd[1]
  if(nrow(sub.table)>1){ # if there are other matches
    for(kk in 2:nrow(sub.table)){
      set.pos <- union(set.pos, sub.table$QStart[kk]:sub.table$QEnd[kk])
    }
  }
  return(length(set.pos))
}

remove.subset <- function(sub.table){ # remove the blast hits that are a subset of longer hits
  idx.to.remove <- c()
  for(ii in 1:(nrow(sub.table)-1)){
    for(jj in (ii+1):nrow(sub.table)){
      if(all(sub.table$QStart[ii]:sub.table$QEnd[ii] %in% sub.table$QStart[jj]:sub.table$QEnd[jj])){
        idx.to.remove <- c(idx.to.remove, ii)
      } else {
        if(all(sub.table$QStart[jj]:sub.table$QEnd[jj] %in% sub.table$QStart[ii]:sub.table$QEnd[ii]))idx.to.remove <- c(idx.to.remove, jj)
      }      
    }
  }
  return(sub.table[setdiff(1:nrow(sub.table), idx.to.remove),])
  # TO DO issue a warning if hits are in opposite directions
}

#### functions to calculate the complement of a sequence and the reverse complement ####

complement <- function(nuc){
  if(nuc == "-")return("-")
  if(nuc == "a")return("t")
  if(nuc == "c")return("g")
  if(nuc == "g")return("c")
  if(nuc == "t")return("a")
  stop("nucleotide not supported")
}
reverse.complement <- function(seq){ # seq in vector format (one entry for each nucleotide)
  return(sapply(rev(seq), complement))
}


get.bee.id <- function(label){
  return(substr(label,1,7))
}

if(!"seqinr" %in% rownames(installed.packages()))install.packages("seqinr", repos="http://cran.rstudio.com/")
library(seqinr)

# read the contigs and the blast matches
all.contigs <- read.fasta(contigs.filename, as.string = F)
all.names <- names(all.contigs)

blast.matches <- read.csv(blast.filename, header=F)
# columns of the blast output are:
# qacc means Query accesion
# sacc means Subject accession
# sseqid means Subject Seq-id
# evalue means Expect value
# pident means Percentage of identical matches
# qstart means Start of alignment in query
# qend means End of alignment in query
# sstart means Start of alignment in subject
# send means End of alignment in subject

stopifnot(ncol(blast.matches) == 9) # check 9 columns

blast.matches <- blast.matches[,-3] # remove third column
colnames(blast.matches)<-c("Contig_name","Match","Eval","PropI","QStart","QEnd","MStart","MEnd")
cat(nrow(blast.matches), " blast matches") # 2796 matches of 24891 contigs recovered
cat(" corresponding to ", length(unique(blast.matches$Contig_name)), " unique contigs\n")
all.id <- unique(sapply(blast.matches$Contig_name, get.bee.id))
cat(length(unique(all.id)), " patients\n") # 782 patients


cleaned.contigs <- sort.blast.output(blast.matches, all.contigs)
cleaned.contigs.seq <- cleaned.contigs[[1]]
cleaned.contigs.names <- cleaned.contigs[[2]]
id.post.cleaning <- sapply(cleaned.contigs.names, get.bee.id)
cat(length(unique(id.post.cleaning)), " patients post contig cleaning\n") # still 782 patients

#### write the results in files ####
cat("writing cleaned alignment in file...\n")
order(cleaned.contigs.names) -> reorder
write.fasta(sequences=cleaned.contigs.seq[reorder],
            names = cleaned.contigs.names[reorder],
            file.out=alignment.out.name,nbchar=80)


