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
    help="Rainstorm output (points to plot) from a single chromosome",default="rainstorm_noncoding_k_4_mean_3.tsv"
    );
parser$add_argument(
    "--xstart", "-s",
    help="limit plot to points after this position",default=-1
);
parser$add_argument(
    "--xend", "-e",
    help="limit plot to points before this position",default=-1
        );


args = parser$parse_args();

infile = args$rainstorm_points;



xmin = as.integer(args$xstart);
xmax=as.integer(args$xend)

if(args$xend<0){
xmax = NA
}else{
xmax = as.integer(args$xend);
}
if(args$xstart<0){
xmin = NA
}else{
xmin=as.integer(args$xstart)
}
print(paste("using","xmin",xmin,"xmax",xmax));
points=read.csv(infile,sep="\t",stringsAsFactors=F);
#this option shows the points above 0, which are largely uninformative
#ggplot(points,aes(x=position,y=mutrate,colour=patient)) + geom_point(alpha=0.2) + theme(legend.position="none")
if(is.na(xmin) && is.na(xmax)){
    ggplot(points,aes(x=position,y=mutrate,colour=patient)) + geom_point(alpha=0.2) + theme_classic() + theme(legend.position="none") + ylim(NA,0) 
}else{
    ggplot(points,aes(x=position,y=mutrate,colour=patient)) + geom_point(alpha=0.2) + theme_classic() + theme(legend.position="none") + ylim(NA,0) + xlim(xmin,xmax)
    }


ggsave(file=paste(infile,".pdf",sep=""),width=7,height=4)
