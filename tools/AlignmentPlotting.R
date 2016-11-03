library(ggplot2)
library(reshape)
library(grid)
library(gridExtra)
library(gtable)
library(scales)
library(argparse)

arg_parser = ArgumentParser()

arg_parser$add_argument("covFile", metavar="coverageFileName")  
arg_parser$add_argument("colFile", metavar="coloursFileName")  
arg_parser$add_argument("genFile", metavar="genesFileName")  
arg_parser$add_argument("outFile", metavar="outputileName")  

args <- arg_parser$parse_args()

cov.file.name <- args$covFile
col.file.name <- args$colFile
gen.file.name <- args$genFile
out.file.name <- args$outFile

AlignPlots <- function(...) {
  LegendWidth <- function(x) x$grobs[[8]]$grobs[[1]]$widths[[4]]
  
  plots.grobs <- lapply(list(...), ggplotGrob)
  
  max.widths <- do.call(unit.pmax, lapply(plots.grobs, "[[", "widths"))
  plots.grobs.eq.widths <- lapply(plots.grobs, function(x) {
    x$widths <- max.widths
    x
  })
  
  legends.widths <- lapply(plots.grobs, LegendWidth)
  max.legends.width <- do.call(max, legends.widths)
  plots.grobs.eq.widths.aligned <- lapply(plots.grobs.eq.widths, function(x) {
    if (is.gtable(x$grobs[[8]])) {
      x$grobs[[8]] <- gtable_add_cols(x$grobs[[8]],
                                      unit(abs(diff(c(LegendWidth(x),
                                                      max.legends.width))),
                                           "mm"))
    }
    x
  })
  
  plots.grobs.eq.widths.aligned
}

process.string <- function(string){
  anything.exists <- F
  
  gs <- gregexpr("g*", string)
  lengths <- attr(gs[[1]], "match.length")
  nonzero.lengths <- which(lengths!=0)
  if(length(nonzero.lengths)>0){
    starts <- gs[[1]][nonzero.lengths]
    display.starts <- starts - 0.5
    ends <- starts + lengths[nonzero.lengths] - 1
    display.ends <- ends + 0.5
    out.g <- data.frame(start = starts, end = ends, display.start = display.starts, display.end = display.ends, char = "g")
    
    anything.exists <- T
    out <- out.g
  }
  
  ds <- gregexpr("d*", string)
  lengths <- attr(ds[[1]], "match.length")
  nonzero.lengths <- which(lengths!=0)
  if(length(nonzero.lengths)>0){
    starts <- ds[[1]][nonzero.lengths]
    display.starts <- starts - 0.5
    ends <- starts + lengths[nonzero.lengths] - 1
    display.ends <- ends + 0.5
    out.d <- data.frame(start = starts, end = ends, display.start = display.starts, display.end = display.ends, char = "d")
    
    if(anything.exists){
      out <- rbind(out, out.d)
    } else {
      anything.exists <- T
      out <- out.d
    }
  }
  
  bs <- gregexpr("b*", string)
  lengths <- attr(bs[[1]], "match.length")
  nonzero.lengths <- which(lengths!=0)
  if(length(nonzero.lengths)>0){
    starts <- bs[[1]][nonzero.lengths]
    display.starts <- starts - 0.5
    ends <- starts + lengths[nonzero.lengths] - 1
    display.ends <- ends + 0.5
    out.b <- data.frame(start = starts, end = ends, display.start = display.starts, display.end = display.ends, char = "b")
    
    if(anything.exists){
      out <- rbind(out, out.b)
    } else {
      anything.exists <- T
      out <- out.b
    }
  }
  
  return(out)
}


cov.data <- read.table(cov.file.name, sep=",", header=T, stringsAsFactors = F)
cov.data <- cov.data[,c(1,4,5)]
comp.factors <- c(colnames(cov.data)[c(3,2)])
cov.data <- melt(cov.data, id = 1)
cov.data$variable <- factor(cov.data$variable, comp.factors)

pos.data <- read.table(col.file.name, sep=",", header=F, stringsAsFactors = F)
colnames(pos.data) <- c("Name", "Sequence")
pos.data$Name <- factor(pos.data$Name, rev(pos.data$Name))

line.list <- lapply(pos.data$Sequence, process.string)

first <- T

for(sequence.no in seq(1, nrow(pos.data))){
  temp <- line.list[[sequence.no]]
  temp$phantom.pos <- sequence.no
  temp$sequence.name <- pos.data$Name[sequence.no]
  if(first){
    line.df <- temp
    first <- F
  } else {
    line.df <- rbind(line.df, temp)
  }
  
}

ann.data <- read.table(gen.file.name, sep=",", header=T, stringsAsFactors = F)
ann.data$display.start <- ann.data$start-0.5
ann.data$display.end <- ann.data$end+0.5
ann.data <- ann.data[order(ann.data$gene),]

ann.data.main <- ann.data[which(ann.data[,5]=="no"),]
ann.data.extra <- ann.data[which(ann.data[,5]=="yes"),]

graph.cov <- ggplot(data=cov.data)

graph.cov <- graph.cov + 
  geom_line(aes(x=Alignment.position, y=value, colour=variable)) +
  ylab("Coverage") +
  xlab("Alignment position") +
  scale_x_continuous(limits=c(min(cov.data$Alignment.position)-1, max(cov.data$Alignment.position)+1), expand=c(0,100)) +
  scale_y_log10(limits = c(1,NA), labels = comma) +
  scale_colour_manual(labels=c("HXB2", "shiver reference"), values=c("red", "blue")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm")) 

graph.ann <- ggplot()

graph.ann <- graph.ann + 
  geom_segment(data = ann.data.main, aes(y=Reading.frame, yend = Reading.frame, x = display.start, xend=display.end, colour=gene), size=5) +
  geom_text(data = ann.data.main, aes(label=gene, y=Reading.frame, x=display.start), hjust="left", nudge_x=10) +
  geom_segment(data = ann.data.extra, aes(y=Reading.frame, yend = Reading.frame, x = display.start, xend=display.end), size=5, alpha=0.5, colour="grey10") +
  geom_text(data = ann.data.extra, aes(label=gene, y=Reading.frame, x=display.start), hjust="left", nudge_x=10, colour="white") +
  theme_bw() +
  scale_x_continuous(limits=c(min(cov.data$Alignment.position)-1, max(cov.data$Alignment.position)+1), expand=c(0,100)) +
  scale_y_continuous(limits=c(0, max(ann.data$Reading.frame)+1)) +
  scale_color_discrete(h = c(0, 720) + 15, direction = -1) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        plot.background=element_blank(),
        plot.margin=unit(c(0.5,0.5,0,0.5), "cm"),
        panel.margin=unit(c(0.5,0.5,0,0.5), "cm"))

graph.pos <- ggplot(data = line.df)

graph.pos <- graph.pos +
  geom_segment(aes(y=sequence.name, yend = sequence.name, x=display.start, xend = display.end, colour=char, size=char)) +
  theme_bw() +
  scale_size_manual(values=c(5,0.5,5)) +
  scale_color_manual(values=c("grey85", "black", "black")) +
  scale_x_continuous(limits=c(min(cov.data$Alignment.position)-1, max(cov.data$Alignment.position)+1), expand=c(0,100)) +
  theme(legend.position="none", 
        axis.ticks.y=element_blank(), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(), 
        plot.margin=unit(c(0,0.5,0,0.5), "cm"),
        panel.margin=unit(c(0,0.5,0,0.5), "cm"))


plots1 <- AlignPlots(graph.ann, graph.pos, graph.cov)
plots1$ncol <- 1
plots1$heights <- unit(c(2,3,6), "null")

pdf(file=out.file.name, width=25, height=10)
do.call(grid.arrange, plots1)
dev.off()