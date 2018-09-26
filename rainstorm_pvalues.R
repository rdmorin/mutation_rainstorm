##################################################################
####### R code for calculating empirical P-values for each wavelet peak
####### Aixiang Jiang, 20180705
###################################################################

### libraries
suppressMessages(library(argparse))
suppressMessages(library(maftools))
suppressMessages(library(data.table))
suppressMessages(library(bedr))


#### input parameters
parser <- ArgumentParser(description="wavelet peak empirical testing argument");

parser$add_argument(
  "--input_maf", "--m",
  help="MAF file containing mutation calls from many patient genomes"
);

parser$add_argument(
  "--peak_summary", "-peak", help="Wavelet peak summary file" 
);


parser$add_argument(
  "--gap_file", "-gap", help="Mutation gap position ranges for all Chromosomes" 
);

parser$add_argument(
  "--output_base_name","--o",help="specify a base file name prefix for all outputs"
);

parser$add_argument(
  "--perm_N","--N",help="specify number of permutations", default = 100000
);

args = parser$parse_args();

maf.file=args$input_maf
peakSumFile = args$peak_summary
gapFile = args$gap_file
outname = args$output_base_name
perm_n = as.integer(args$perm_N)

#########################################
#### read in maf file
maf.full = fread(maf.file, verbose = FALSE, key = c("Chromosome", "Start_Position", "End_Position"))

#### read in peak summary data
peaks = read.delim(peakSumFile, header=T, row.names = 1, stringsAsFactors = F)

##### read in gap data
gaps = read.delim(gapFile, header=T, stringsAsFactors = F)
tgaps = data.table(gaps, key = c("chrom", "start", "end"))


######################################

getMafRange = function(mafdat){
  ### I need chr, and for each chr, get start and end positions
  chrs = paste("chr", c(1:22,"X"), sep="")
  chrout = t(sapply(chrs, FUN = function(chrin){
    schr = subset(mafdat, mafdat$Chromosome == chrin)
    c(min(schr$Start_Position), max(schr$End_Position))
  }))
  return(cbind(chrs, chrout))
}

chrinfo = getMafRange(mafdat = maf.full)
colnames(chrinfo) = c("chr","start","end")
chrinfo = data.frame(chrinfo, stringsAsFactors = F)
chrinfo$start = as.integer(chrinfo$start)
chrinfo$end = as.integer(chrinfo$end)


#########################
tchrinfo = data.table(chrinfo, key = c("chr","start","end"))
chrRmGap = foverlaps(tgaps, tchrinfo, by.x =c("chrom", "start", "end"), by.y = c("chr","start","end"))
#### this is not what I want exactly, remove Y
chrRmGap = chrRmGap[-c(25:26),]
chrRmGap = chrRmGap[, c(1,2,4,5,3)]
chrRmGap$i.start = as.integer(ifelse(chrRmGap$i.start >= chrRmGap$start, chrRmGap$i.start, chrRmGap$start+1))
tmp21 = subset(chrRmGap, chrRmGap$chrom == "chr21")
tmp21$end[1] = as.integer(mean(c(tmp21$i.end + tmp21$i.start)))
tmp21$start[2] = tmp21$end[1]+1

tmp = which(chrRmGap$chrom == "chr21")
tmp21$start = as.integer(tmp21$start)
chrRmGap[tmp,] = tmp21

chlength = apply(chrRmGap[,-1], 1, FUN = function(xx){
  sizes = length(c(xx[1] : (xx[2]-1), (xx[3]+1):xx[4]))
})

totlength = sum(as.numeric(chlength))
props = chlength/totlength
chrall = cbind(chrRmGap, chlength, props)
colnames(chrall)[c(3:4,6:7)] = c("gapStart","gapEnd","validLength","propToGenome")

##### a function to find p value for each peak
getEmpiricalP = function(apeak, permN=100000, mafDatTable = maf.full, mafranges = chrall){
 
  peaksize = as.integer(apeak[3]) - as.integer(apeak[2]) +1 
  
  mchoise = dim(mafranges)[1]
  
  chrsim = sample(1:mchoise, size = permN, replace = T)
  mafrangeInd = data.frame(table(chrsim))
  mafrangeInd$chrsim = as.integer(mafrangeInd$chrsim)
  
  chrout = apply(mafrangeInd,1,FUN = function(rangeIndSize){
    rangeInd = rangeIndSize[1]
    rangeFreq = rangeIndSize[2]
    simRange = c(as.integer(mafranges[rangeInd,"start"]):as.integer(mafranges[rangeInd,"gapStart"]-1), 
                 as.integer(mafranges[rangeInd,"gapEnd"]+1):as.integer(mafranges[rangeInd,"end"]))
    asim = sample(simRange, size = rangeFreq)

    asim1 = asim - as.integer(peaksize/2)
    asim2 = asim + as.integer(peaksize/2) -1
    
    chrInfo = as.character(rep(mafranges[rangeInd,1], rangeFreq))
    
    return(cbind(chrInfo, asim1, asim2))
  })
  
  chrout = do.call(rbind, chrout)
  colnames(chrout) = c("chrom","start", "end")

  simIn = data.frame(chrout, stringsAsFactors = F)
  simIn$start = as.integer(simIn$start)
  simIn$end = as.integer(simIn$end)
  tsimIn = data.table(simIn, key = c("chrom","start", "end"))
  tmatchs = foverlaps(mafDatTable, tsimIn, by.x = c("Chromosome","Start_Position", "End_Position") , by.y = c("chrom","start", "end"))
  tmatchs = tmatchs[!is.na(tmatchs$start),]
  muts = table(tmatchs$start)
  mutrates = muts*1000/(peaksize)
  mm = dim(tsimIn)[1] - length(muts)
  mutrates = c(mutrates, rep(0,mm))
  peakmut = as.numeric(apeak[4])
  pvalue = length(which(mutrates >= peakmut))/permN
  
  pdffile = paste(toString(apeak),"_Pvalue",pvalue,"_",Sys.Date(),".pdf",sep="")
  pdffile = gsub(",","_", pdffile)
  pdffile = gsub(" ","", pdffile)
  pdffile = paste(outname, pdffile, sep="")
  pdf(pdffile)
  hist(mutrates, breaks = 100)
  dev.off()
  
  return(c(peakmut, permN,pvalue))
}

#### get four columns from peaks for the testing: chr LeftPosition RightPosition  mutPerKb
#### in order to work for outputs from both wavelet R scripts, do the following
t1 = grep("chr",colnames(peaks))
t2 = grep("eftPosition",colnames(peaks))
t3 = grep("ightPosition",colnames(peaks))
t4 = which(colnames(peaks)=="mutPerKb")
peaksIn = peaks[,c(t1,t2,t3,t4)]

getpermP = t(apply(peaksIn, 1, getEmpiricalP, permN=perm_n, mafDatTable = maf.full, mafranges = chrall))
colnames(getpermP) = c(c("peak_mutrate", "permutationN", "empirical_Pvalue"))
outall = cbind(peaks, getpermP[rownames(peaks),])

write.table(outall, paste(outname,"_peak_empirical_test_",Sys.Date(),".tsv",sep=""), sep="\t", quote = F)
