##########################################
##### Rainstorm calculation - cohort-wide variation of rainfall plots based on a MAF file from many cancer genomes
##### Authors: Ryan Morin and Aixiang Jiang, 2017
########################################

library(argparse)

require(ggplot2)
library(data.table)



parser <- ArgumentParser(description="plot a genomic region with mutations and encode data");


parser$add_argument(
    "--file_base_name","-bn",
    help="characters before chr# or # in the input file names (required when plotting values from more than one chromosome)"
    );

parser$add_argument(
    "rainstorm_files",
    help="Full path to one or more rainstorm output files", nargs='+', type="character"
    );

parser$add_argument(
           "--chromosome_name","-c",
           help="optionally specify a single chromosome to plot",
           default="all");

parser$add_argument(
    "--xstart", "-s",
    help="limit plot to points after this position",default=-1
);
parser$add_argument(
    "--xend", "-e",
    help="limit plot to points before this position",default=-1
        );
parser$add_argument(
           "--plot_style","-P",
           help="use a circular or linear axis for your chromosomal coordinate?",
           default="linear"
           );

args = parser$parse_args();

infile = args$rainstorm_points;
plotstyle = args$plot_style;
chromname = args$chromosome_name;
chromfile = args$rainstorm_chr;
allfiles = args$rainstorm_files;
xmin = as.integer(args$xstart);
xmax=as.integer(args$xend)
separator = args$file_base_name;



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

    

load_raw_points <-  function(files=NULL,separator=NULL){
    points.list = list()
    for(f in files){
        trimmed = gsub(".tsv","",f)
        chrom = gsub(separator,"",trimmed)
        print(paste(f,chrom));
        points=read.csv(f,sep="\t",stringsAsFactors=F);
        points$chromosome = chrom;
        points$relative_pos = points$position / max(points$position);
        points.list[[chrom]]=points
    }
    points.all = do.call("rbind",points.list);
}

plot_single_chrom <- function(style="linear",chrom=NULL,base=NULL){
    if(style == "circular"){
        ggplot(points,aes(x=position,y=mutrate,colour=patient)) + geom_point(alpha=0.2,size=0.2) + theme_classic() + theme(legend.position="none") + ylim(NA,0) + coord_polar(theta ="x");
     ggsave(file=paste(base,chrom,style,"_single",".png",sep=""),width=9,height=9)
        ggsave(file=paste(base,chrom,style,"_single",".pdf",sep=""),width=9,height=9)
        
    } else{
        
        if(is.na(xmin) && is.na(xmax)){
            ggplot(points,aes(x=position,y=mutrate,colour=patient)) + geom_point(alpha=0.2,size=0.2) + theme_classic() + theme(legend.position="none") + ylim(NA,0) 
        }else{
            ggplot(points,aes(x=position,y=mutrate,colour=patient)) + geom_point(alpha=0.2,size=0.2) + theme_classic() + theme(legend.position="none") + ylim(NA,0) + xlim(xmin,xmax)
        }
        ggsave(file=paste(base,chrom,style,"_single",".png",sep=""),width=16,height=3)
        ggsave(file=paste(base,chrom,style,"_single",".pdf",sep=""),width=16,height=3)
    }
    

}

plot_all_chrom <-function(style="linear",basename=separator){
    if(style == "circular"){
        ggplot(points,aes(x=relative_pos,y=mutrate,colour=patient)) + geom_point(alpha=0.1,size=0.1) + theme_classic() + theme(legend.position="none") + ylim(NA,0) +
        coord_polar(theta="x") + facet_wrap(~chromosome,ncol=5); #plot in a 5x5 grid 
        ggsave(file=paste(basename,"_",style,"_all_chr.pdf",sep=""),width=18,height=18)
    }else{
        ggplot(points) + geom_point(aes(colour=patient,x=position,y=mutrate),alpha=0.1,size=0.1) + facet_wrap(~chromosome,scale="free_x",ncol=2) +
            theme_classic() + theme(legend.position = 'none') + ylim(NA,0)
        print(paste("saving as",basename,"allchr_linear.pdf"))
        ggsave(file=paste(basename,"allchr_linear.pdf",sep=""),width=14,height=22)
        ggsave(file=paste(basename,"allchr_linear.png",sep=""),width=14,height=22)
    }
}


if(length(allfiles)>0){
    points = load_raw_points(allfiles,separator);
}else{
    points=read.csv(infile,sep="\t",stringsAsFactors=F);
    }
print(paste("chrom:",chromname))
if(chromname == "all"){
    plot_all_chrom(plotstyle);
}else{
    plot_single_chrom(plotstyle,chrom=chromname,separator);
}
