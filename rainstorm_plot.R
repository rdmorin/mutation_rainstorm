##########################################
##### Rainstorm calculation - cohort-wide variation of rainfall plots based on a MAF file from many cancer genomes
##### Authors: Ryan Morin and Aixiang Jiang, 2017
########################################

library(argparse)

require(ggplot2)
library(data.table)



parser <- ArgumentParser(description="plot a genomic region with mutations and encode data");

parser$add_argument(
    "--rainstorm_points", "-r",
    help="points from a single chromosome",default="rainstorm_noncoding_k_4_mean_3.tsv"
    );
parser$add_argument(
    "--xstart", "-s",
    help="limit plot to points after this position"
);
parser$add_argument(
    "--xend", "-e",
    help="limit plot to points before this position"
        );


args = parser$parse_args();

infile = args$rainstorm_points;
xmin = NA;
xmax = NA;
if(!is.na(args$xstart)){
    xmin = as.integer(args$xstart);
}
if(!is.na(args$xend)){
    xmax = as.integer(args$xend);
}
points=read.csv(infile,sep="\t",stringsAsFactors=F);
#this option shows the points above 0, which are largely uninformative
#ggplot(points,aes(x=position,y=mutrate,colour=patient)) + geom_point(alpha=0.2) + theme(legend.position="none")
ggplot(points,aes(x=position,y=mutrate,colour=patient)) + geom_point(alpha=0.2) + theme_classic() + theme(legend.position="none") + ylim(NA,0) + xlim(xmin,xmax)


ggsave(file=paste(infile,".pdf",sep=""),width=7,height=4)
