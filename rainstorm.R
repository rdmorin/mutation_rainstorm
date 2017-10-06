##########################################
##### Rainstorm calculation - cohort-wide variation of rainfall plots based on a MAF file from many cancer genomes
##### Authors: Ryan Morin and Aixiang Jiang, 2017
########################################

library(argparse)
require(GenomicRanges)
require(ggplot2)
library(maftools)
library(data.table)
library(parallel)


parser <- ArgumentParser(description="plot a genomic region with mutations and encode data");

parser$add_argument(
    "--input_maf", "-m",
    help="MAF file containing mutation calls from many patient genomes"
    );
parser$add_argument(
           "--cpu_num","-c",help="set to number of CPUs you would like to use to perform calculation in parallel (consumes lots of RAM)",default=1);
parser$add_argument(
    "--genome_fai", "-g",
    help="provide the corresponding fasta index for the genome you used. must match the chromosome naming style used in your MAF!", default="hg19.ensembl.fa.fai"
);

parser$add_argument("--plot","-p",help="ploduce rainstorm plot for each chromosome",default=TRUE);
parser$add_argument(
"--max_mut","-M",help="genomes skipped if their total mutation load exceeds this value",default=50000
);
parser$add_argument(
           "--off_by","-k",help="take mean of the distance to the k closest mutations to determine rainstorm distance value", default=4)

parser$add_argument(
    "--calc_background","-b",help="if you have done this once for a cohort, you can reload the result in future runs by setting this to 0",default=1);

args = parser$parse_args();

genome.fai = args$genome_fai
off.by = as.integer(args$off_by);
cpu.num=as.integer(args$cpu_num);
calc.background = as.integer(args$calc_background);

mutcount.max = as.integer(args$max_mut);

if(!is.null(genome.fai)){
    genomedetails = read.table(genome.fai,sep="\t")
    print(genomedetails);
    chrlengths  = genomedetails[,2];
    names(chrlengths) = genomedetails[,1];

}else{
    print("error, missing genome");
    q();
}
maf.file = args$input_maf

maf.full = read.maf(maf.file,useAll = T, removeSilent = F)

#get IDs of cases passing the max mutation criteria
use.cases = as.character(maf.full@variants.per.sample[Variants<mutcount.max,Tumor_Sample_Barcode])
                                       


#latest version of the function to correct local mutation rate. This uses a loess model that fits a smoothed curve to the mutation rate across the chromosome
correctLocalMutrate2<-function(chrom,positions,distval,model,logged_mutrate){
  predrate =  predict(model,newdata=data.frame(starts=positions))
  adjusted = log(as.numeric(distval)) + predrate + logged_mutrate
  return(adjusted)
}

#function used to obtain a value for each mutation that is later scaled for local mutation rate. This function just compares pairs of genomes and is called by another function when performing a one-to-all comparison
getMutDists <-function(pos1,pos2,id1='G1',getmin=FALSE){
  if(length(pos1)==0 || length(pos2)==0){
    return(rep("NA",length(pos1)))
  }
  id2 = "G2"
  these = data.frame(positions=c(pos1,pos2),genomes=factor(c(rep(id1,length(pos1)),rep(id2,length(pos2)))))
  
  these = these[order(these[,1]),]  #sort on position
  if(getmin==FALSE){
    these[,"dist"] = diff(c(these[,1],NA))
    
    these[these[,"genomes"]==id1,"distself"]=diff(c(these[these[,"genomes"]==id1,"positions"],NA))
    these[these[,"genomes"]==id2,"distself"]=diff(c(these[these[,"genomes"]==id2,"positions"],NA))
    
    these[,"keepdist"]=these[,"dist"]
    mask= which(these[,"dist"]==these[,"distself"])
    these[mask,"keepdist"] = NA
  }else{
    these[,"dist"] = diff(c(these[,1],NA))
    these[,"distrev"] = diff(c(NA,these[,1]))
    these[these[,"genomes"]==id1,"distself"]=diff(c(these[these[,"genomes"]==id1,"positions"],NA))
    these[these[,"genomes"]==id2,"distself"]=diff(c(these[these[,"genomes"]==id2,"positions"],NA))
    
    these[these[,"genomes"]==id1,"distselfrev"]=diff(c(NA,these[these[,"genomes"]==id1,"positions"]))
    these[these[,"genomes"]==id2,"distselfrev"]=diff(c(NA,these[these[,"genomes"]==id2,"positions"]))
    
    these[,"keepdist"]=  as.numeric(apply(these,1,function(x){min(x["distrev"],x["dist"],na.rm=T)}))
    mask= which(these[,"keepdist"]==these[,"distself"] | these[,"keepdist"]==these[,"distselfrev"])  #this logic would miss hot spots and needs to be fixed. at a hot spot, two mutations will have the same distance to their neighbour
    these[mask,"keepdist"] = NA
    
    #mins = apply(cbind(d,d1),1,min) 
  }
  return(these[these[,"genomes"]==id1,"keepdist"]) #distance to nearest SNV in the other genome or NA if the nearest SNV is in this genome
  
}

#calls the getMutDists function on all cases for a single index case (id)
getMinDistByGenome<-function(maf,id,chromosome,use.cases,start,end,offby=2,getmin=FALSE,usemean=FALSE){
  #extract mutations in region for this genome and compute the N-closest minimum distance to each variant among all genomes (default N, 2). Self is ignored, nearest genome is ignored.
  ### back to noncoding only
  thesemut = as.data.frame(maf@data[Tumor_Sample_Barcode==id & Variant_Classification %in% noncoding & Chromosome == chromosome & Start_Position > start & Start_Position < end,Start_Position])[,1]

    thesemut = thesemut[order(thesemut)]
  distmat = matrix(nrow=length(thesemut),ncol=length(use.cases)-1)
  i=1
  for(mycase in use.cases){
    if(mycase == id){
      next; 
    }
    thosemut = as.data.frame(maf@data[Tumor_Sample_Barcode==mycase & Chromosome == chromosome &  Start_Position > start & Start_Position < end,Start_Position])[,1]
    thosemut = thosemut[order(thosemut)]
    dists =  getMutDists(thesemut,thosemut,getmin=getmin) #getmin option (not used currently) would force the distance in both directions to be considered and keep min()
    distmat[,i]=dists
    i=i+1
  }
  if(usemean){
    keepdist = apply(distmat,1,function(x){mean(x[order(x)][c(1:offby)])})
  }else{
      keepdist = apply(distmat,1,function(x){x[order(x)][offby]})
      #original approach is to just return the kth value instead of mean from 1:k
  }
  return(data.frame(position=thesemut,mindist = keepdist,stringsAsFactors = F))
}
 
noncoding = c("3'Flank","IGR","Intron","3'UTR","5'Flank","5'UTR","Targeted_Region","RNA")



############################################################################################################
### start from here: for binned data, which should not be re-run for the same disease data type ##############

#goodchrom =gsub("chr","",names(goodchrom.len))
goodchrom =names(chrlengths) #may want to add an option where user can specify the list of chromosomes to include
snvs = GRanges(seqlengths=chrlengths,
                  ### back to noncoding
                      seqnames=maf.full@data[ Variant_Classification %in% noncoding & Chromosome %in% goodchrom &
                                              Tumor_Sample_Barcode %in% use.cases,Chromosome],
                      IRanges(maf.full@data[  Variant_Classification %in% noncoding &
                                              Chromosome %in% goodchrom & Tumor_Sample_Barcode %in% use.cases,Start_Position],width=1),
                      names=maf.full@data[Variant_Classification %in% noncoding & Chromosome %in% goodchrom &
                                          Tumor_Sample_Barcode %in% use.cases,Tumor_Sample_Barcode])


if(calc.background){
#100kb bins, make size an option?
bins.chr = tileGenome(seqlengths=chrlengths,tilewidth = 100000)

bincounts.all = c()
binstarts.all=c()
bincounts.chrom=c()
binstops.all = c()
bin_length = 100000
for(chrom in goodchrom){

  print(chrom)
  patient = use.cases[1]
  chr = chrom

  cvg <- coverage(snvs[snvs$names==patient,])
  pat.tot=length(snvs[snvs$names==patient,]) #each patient's values are reduced to reflect the genome-wide mutation rate
  r <- runsum(cvg, 1)
  tile=binnedAverage(unlist(bins.chr),r,"binned_score")
  ntile =length(tile[seqnames(tile)==chr,]$binned_score )

  npat=length(use.cases)
  testmat = matrix(nrow=ntile,ncol=npat)
  for(num in c(1:npat)){

    patient  = use.cases[num]
    cvg <- coverage(snvs[snvs$names==patient,])
    pat.tot=length(snvs[snvs$names==patient,])
    r <- runsum(cvg, 1)

    tile=binnedAverage(unlist(bins.chr),r,"binned_score")

    a=tile[seqnames(tile)==chr,]$binned_score
    testmat[,num]=a[c(1:ntile)]
  }
  means = apply(testmat,1,mean)
  means = means * 0.000000001
  logmeans = log(means*bin_length) #note, natural log scale (harmonize to log10 as per the distance/rainfall approach?)
  bincounts.all=c(bincounts.all,logmeans)
  binstarts.all=c(binstarts.all,start(tile[seqnames(tile)==chr,]))
  binstops.all=c(binstops.all,end(tile[seqnames(tile)==chr,]))

  bincounts.chrom=c(bincounts.chrom,rep(chrom,length(logmeans)))

}
print("done calculating background correction")
#load this next file in as a data frame to skip the steps leading up to this line
alldf = data.frame(chrom=bincounts.chrom,starts=binstarts.all,ends=binstops.all,counts=bincounts.all)

write.table(alldf,file="./background_100k_binned_noncodingsnv_density.tsv",sep="\t",quote=F)
model= loess(counts~starts,data=alldf[alldf[,1]== 3 & alldf[,"counts"]!=-Inf,],span=0.005,surface='direct')
chrdf =alldf[alldf[,1]==3,]
lopoints=predict(model,newdata=chrdf)

chrdf$predicted=lopoints
ggplot(chrdf,aes(x=starts,y=counts)) + geom_point(alpha=0.4,colour='orange') + geom_line(aes(x=starts,y=predicted),colour='red') + ylim(-18,-23); #pllot out chrom3
ggsave(file="chr3_background.pdf",width=7,height=4)

}else{
alldf = read.csv("./background_100k_binned_noncodingsnv_density.tsv",sep="\t",stringsAsFactors=F);
    }

######### end of binned data part, which should not be re-run for the same disease data type ######################
#####################################################################################################################


##### move functions here, before calling them

runByCaseSmooth<-function(case,maf.full,background.model,use.cases,chrom,start,end,offby=3,getmin=FALSE){
  stored.all = list()
  use.mean=TRUE
  these = getMinDistByGenome(maf.full,case,chrom,use.cases,start,end,offby=offby,getmin=getmin,usemean=use.mean)
  if(dim(these)[1]==0){
    print(paste("skip due to lack of mutations on chromosome",case)) 
    
    return(stored.all)
  }
  densitytable = table(is.na(these))
  if(length(densitytable)>1){

    if(densitytable[2]/densitytable[1] > nathresh){
      print(paste("skip due to high NA count",case))

      return(stored.all)
    }
  }
  
  tot = length(these[,1])
  genometot = maf.full@variants.per.sample[Tumor_Sample_Barcode==case,Variants]
  ltot = log(genometot/280000000)
  #print(paste("shifting by",genometot,ltot))
  napos = is.na(as.numeric(these[,2]))
  nonatot = length(these[!napos,1])
  these.keep = these[!napos,]
  #should we get rid of all NA values first? Seems reasonable since they're being counted here in the denominator? Though they are in fact mutations, so maybe not...
  scaled = log(as.numeric(these.keep[,"mindist"])+1) + log(genometot)-log(280000000) #add one to get rid of the annoying -Inf issue. These are definitely things that need to be retained.
  scaled = scaled - median(scaled)
  localadj = correctLocalMutrate2(chrom,these.keep[,1],as.numeric(these.keep[,"mindist"])+1,background.model,ltot)
  
  scaled.localadj=localadj - median(localadj,na.rm=T)
  stored.all[["mutdiff"]]=as.numeric(these.keep[,"mindist"])
  stored.all[["position"]]=these.keep[,1]
  stored.all[["mutrate"]]=scaled.localadj
  stored.all[["mutrate.noadj"]]=scaled
  stored.all[["patient"]]=rep(case,length(these.keep[,1]))
  return(stored.all)
}

plotRainstorm<-function(points,name){
    ggplot(points,aes(x=position,y=mutrate,colour=patient)) + geom_point(alpha=0.2) + theme_classic() + theme(legend.position="none") + ylim(NA,0)
    ggsave(file=name,width=7,height=4)
    }

  
  n = length(use.cases)
  n=n+22
  cols=colors()[c(23:n)]
  names(cols) = use.cases

  nathresh = 0.3 #20% NA values. Need to make this a parameter that users can modify
  
  
  
for(chrom in goodchrom){
    print(paste("running calculation for",chrom));
    start=1
    end = max(maf.full@data[ Chromosome == chrom,Start_Position])
    #chrom = paste("chr",chrom,sep="")
    model= loess(counts~starts,data=alldf[alldf[,1]== chrom & alldf[,"counts"]!=-Inf,],span=0.005,surface='direct')
    if(cpu.num > 1){
    alldat.allpatients = mclapply(use.cases,
      function(x){
      runByCaseSmooth(x,maf.full,model,use.cases,chrom,start,end,offby=off.by)
      },
                                        mc.cores = cpu.num) ### mc.cores is to set number of parallel jobs: ie. CPUs
        }else{
    alldat.allpatients = list();
    for(thiscase in use.cases){
        alldat.allpatients[[thiscase]]=runByCaseSmooth(thiscase,maf.full,model,use.cases,chrom,start,end,offby=off.by)
    }
    }
    #convert to lists with like elements combined and grouped by patient
    patients = c()
    positions = c()
    mutrate = c()
    unadj=c()
    mutdiff=c()
    for(patient in c(1:length(alldat.allpatients))){
        n = length(unadj);
        print(paste(patient,n));
        print("------------");
      unadj = c(unadj,alldat.allpatients[[patient]]$mutrate.noadj)
      patients=c(patients,alldat.allpatients[[patient]]$patient)
      positions=c(positions,alldat.allpatients[[patient]]$position)
      mutrate=c(mutrate,alldat.allpatients[[patient]]$mutrate)
      mutdiff=c(mutdiff,alldat.allpatients[[patient]]$mutdiff)
    }
    #ready points for ggplot rendering
    allcounted = data.frame(mutrate=mutrate,unadj=unadj,position=positions,patient=patients,mutdiff=mutdiff)

    
    
    filen = paste("rainstorm_noncoding_k_",off.by,"_mean_",chrom,".tsv",sep="")
    write.table(allcounted,file=filen,sep="\t",quote=F);
    plotRainstorm(allcounted,gsub(".tsv",".pdf",filen));
    #write out plots here optionally
}
