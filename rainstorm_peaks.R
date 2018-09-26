#######################################################################################
#### R code for peak searching with wavelet, and post processing to modify wavelet peaks
#### Aixiang Jiang, June 2018
########################################################################################


########################### example, Aug 10, 2018 #############################
# Rscript /mnt/thanos_lv/ajiang/AJ2018/RcodeAJ2018/rainWaveTestCodes/waveletPlusPostProcessAug2018.R --output_base_file  /mnt/thanos_lv/ajiang/testwave_FL1 --rainstorm_files  /home/rmorin/rainstorm/FL_smallerbin_rainstorm_k_4_mean_3.tsv --post_process_flag 1 --stringSplit k_4_mean_ --input_maf  /home/rmorin/rainstorm/FL_some_merged.maf
## should change output path if somebody else is used due to folder permission

## check the maf file chr format: /home/rmorin/rainstorm/FL_some_merged.maf
# maf.file = "/home/rmorin/rainstorm/FL_some_merged.maf"
# library(data.table)
# maf.full = fread(maf.file, verbose = FALSE, key = c("Chromosome", "Start_Position", "End_Position"))
# maf.full$Chromosome[1:10]
# > maf.full$Chromosome[1:10]
# [1] "1" "1" "1" "1" "1" "1" "1" "1" "1" "1"
### therefore, have to change the filtering code part for this issue, since bedr needs chr
########################################################

#### there are 11 parameters, and 3 of them do not have default values

#### testing code outside of R in June 2018 
#### under: /mnt/thanos_lv/ajiang/AJ2018/testCombinedWaves
#### Rscript R_wavelet_postProcess_June20CleanVersion.R --rainstorm_files *.tsv --output_base_file testJune20/testJune20 --input_maf /morinlab/projects/NCI_Burkitts/results/tidy/2-paediatric_cohort/strelka.pass.tidy.maf 

suppressMessages(library(argparse))
suppressMessages(library(MassSpecWavelet))
suppressMessages(library(maftools))
suppressMessages(library(data.table))
suppressMessages(library(parallel))
suppressMessages(library(bedr))

########################################################
parser <- ArgumentParser(description="wavelet searching and filtering argument");

parser$add_argument(
  "--rainstorm_files", "-in",  nargs="+",
  help="Rainstorm path and files"  ### this is changed accordingly as well
);

parser$add_argument(
  "--stringSplit", "-split", 
  help="characters before chr# or # in the input file names", default = "_k_4_"
);

parser$add_argument(
  "--input_maf", "-m",
  help="Input maf path and name"
);

parser$add_argument(
  "--output_base_file", "-out", 
  help="Output path and name base"
);

parser$add_argument(
  "--extend_range_1side", "-extend", 
  help="peak extension range for 1 side", default = 10000
);

parser$add_argument(
   "--mut_rate_Per_kb", "-mut_cut", 
   help="mutations per kb threshold", default = 5
);

parser$add_argument(
  "--post_process_flag", "-post", 
  help="Default_1_do__input_0_not_do", default=1
);

parser$add_argument(
  "--distanceToMerge", "-dist", 
  help="maximum distance between regions to be merged", default = 500
);

parser$add_argument(
  "--min_peak_width", "-minWidth", 
  help="min original width", default = 10
);

parser$add_argument(
  "--max_peak_width", "-maxWidth", 
  help="max final peak width", default = 20000
);

parser$add_argument(
  "--patientInRatio", "-petient%", 
  help="min percentage of patients per peak", default = 5 ### actually 5%, need to convert later
);

######################################

args = parser$parse_args();

files=args$rainstorm_files
sepchr = args$stringSplit
outname = args$output_base_file
maf.file = args$input_maf
peakExt1side = as.integer(args$extend_range_1side)
mutRateCut = as.integer(args$mut_rate_Per_kb)
distcut = as.integer(args$distanceToMerge)
postflag = as.integer(args$post_process_flag)
minPeakWidth = as.integer(args$min_peak_width)
maxPeakWidth = as.integer(args$max_peak_width)
noPatCut = as.integer(args$patientInRatio) ### this is in percentage, need to convert later

###############################################################################

#################################################
### get maf file first
maf.full = fread(maf.file, verbose = FALSE, key = c("Chromosome", "Start_Position", "End_Position"))
#### notice that the output is data table, not maf object
####################################################################

set.seed(846161) 
###########################################

#########################################################################
############ functions ################################

getSubMaf = function(awave,mafFileIn){
  names(awave) = c("chromosome", "leftPosition",  "rightPosition")
  awave = data.frame(t(awave), stringsAsFactors = F)
  awave$leftPosition = as.integer(awave$leftPosition)
  awave$rightPosition = as.integer(awave$rightPosition)
  awave = data.table(awave, key =  c("chromosome", "leftPosition",  "rightPosition"))
  
  ### add one line here on July 20
  mafFileIn = data.table(mafFileIn, key = c("Chromosome", "Start_Position", "End_Position"))
  
  subm = foverlaps(x=mafFileIn, y=awave, by.x=c("Chromosome", "Start_Position", "End_Position"), by.y = c("chromosome", "leftPosition",  "rightPosition"))
  subm = subm[!is.na(subm$leftPosition),]
  
  return(subm)
}


getMutPerKb = function(arange, submaf){
  mafsum=subset(submaf, submaf$Start_Position >= as.integer(arange[1]) & submaf$End_Position <= as.integer(arange[2]))
  tn=dim(mafsum)[1]
  mutkb=1000 * tn / ( as.integer(arange[2]) -  as.integer(arange[1])  + 1)
  return(mutkb)
}

getMafSum=function(arange,submaf){  
  ### arange has three items: chr, start and end positions
  mafsum=subset(submaf, submaf$Start_Position >= as.integer(arange[1]) & submaf$End_Position <= as.integer(arange[2]))
  
  tn=dim(mafsum)[1]
  mutkb=1000 * tn / ( as.integer(arange[2]) -  as.integer(arange[1])  + 1)
  vc=mafsum$Variant_Classification
  vc=table(vc)
  vc=names(which.max(vc))
  
  syms=mafsum$Hugo_Symbo
  syms = syms[which(syms != "Unknown" & syms != "")]
  
  syms=table(syms)
  syms=names(which.max(syms))
  
  if(length(syms)==0){syms= "Unknown"}
  
  pts = mafsum$Tumor_Sample_Barcode
  pts = length(unique(pts))
  
  outs = list(c(mutkb,pts),c(vc,syms),mafsum)
  return(outs)
  
}

peakSearch=function(datin, snrperc=0.95, hKmExt =10000, countPerKbcut = 5, chrIn="chr1", mafIn){
  
  ### need chr info as well to map to maf file, and of course, need maf data
  newdat=-datin$mutrate
  peakInfo = peakDetectionCWT(newdat) 
  rm(newdat)
  gc()
  
  majorPeakInfo = peakInfo$majorPeakInfo
  res=cbind(majorPeakInfo$peakSNR,majorPeakInfo$allPeakIndex)
  res=res[which(res[,1]>=0),]
  
  snrCut=quantile(res[,1], snrperc)
  res=res[which(res[,1]>=snrCut),]
  
  #############################################
  
  ridges=peakInfo$ridgeList  
  ridges=ridges[rownames(res)]
  
  ridgeBounds=t(sapply(ridges,FUN=function(x){
    mm=sort(x)
    m1=min(mm)
    m2=max(mm)
    return(c(m1,m2))
  }))
  colnames(ridgeBounds)=c("boundInxMin","boundIndMax")
  
  ridgeBounds=data.frame(cbind(ridgeBounds,res[,1]))
  
  submafIn = subset(mafIn, mafIn$Chromosome == chrIn) 
  submafIn = submafIn[order(submafIn$Start_Position),]  
  
  out=apply(ridgeBounds,1,FUN=function(x){
    resout=NA
    mm=NA
    lm1 = max(c(datin$position[x[1]]-hKmExt, min(datin$position)))
    lm2 = datin$position[x[2]]
    rm1 = datin$position[x[1]]
    rm2 = min(c(datin$position[x[2]]+hKmExt, max(datin$position)))
    
    allsub = getSubMaf(c(chrIn,lm1, rm2 ), mafFileIn = submafIn) 
    omutpkb = getMutPerKb(c(datin$position[x[1]],datin$position[x[2]]),allsub)
    
    if(omutpkb >= countPerKbcut){
      
      ### left side
      l1 = sort(unique(subset(allsub$Start_Position, allsub$Start_Position <= rm1)))
      n1=length(l1)
      l2 = rep(lm2,n1)
      left0 = apply(cbind(l1,l2),1,getMutPerKb, submaf=allsub)
      
      ### right side
      r2 = sort(unique(subset(allsub$End_Position, allsub$End_Position >= lm2)))
      n2 = length(r2)
      r1 = rep(rm1,n2)
      right0 = apply(cbind(r1,r2),1,getMutPerKb, submaf=allsub)
      
      m1=l1
      
      if(length(l1)>1){
        lside = diff(left0)/left0[-1]
        m1 = l1[which.max(lside)+1]
      }
      
      m2=r2
      if(length(r2)>1){
        rside = diff(right0)/right0[-1]
        m2 = r2[which.min(rside)] 
      }
      
      finalMut = getMafSum(c(m1,m2), submaf = allsub)
      ftmp = dim(finalMut[[3]])[1]
      
      if(ftmp>1){
        
        mm=subset(datin, datin$position >= m1 & datin$position <= m2)
        mintmp=which.min(mm$mutrate)
        ntmp=dim(mm)[1]
        
        mm$SNR=rep(x[3], ntmp)
        
        peakpos=mm$position[mintmp]
        peakmin=mm$mutrate[mintmp]
        
        leftpos=min(mm$position)
        rightpos=max(mm$position)
        
        numberP = finalMut[[1]][2]
        
        resout=c(chrIn,peakpos,peakmin,leftpos,rightpos,numberP,mm$SNR[1],finalMut[[1]][1], datin$position[x[1]],datin$position[x[2]], 
                 omutpkb, max(left0),max(right0), finalMut[[2]]) ### there are two items in finalMut[[2]], totally 14 items
        resout=data.frame(t(resout), stringsAsFactors = F)
        colnames(resout) = c("chr","peakMinPosition", "peakMin_mutrate", "LeftPosition", "RightPosition","numberOfPatients","SNR",
                             "mutPerKb","originalLeftPos", "originalRightPos", "original_mutPerKb","max_mutPerKbLeft","max_mutPerKbRight",
                             "mostFreqClassification","mostFreqGeneSymbol")
        resout[,2:13] = apply(resout[,2:13],2,FUN = function(xx) as.numeric(xx))
        
        pn=dim(mm)[1]
        if(pn>1){
          restmp=apply(resout, 2, function(co) rep(co, each = pn))  
        }else{
          restmp=resout
        }
        restmp = data.frame(restmp, stringsAsFactors = F)
        restmp[,2:13] = apply(restmp[,2:13],2,FUN = function(xx) as.numeric(xx))
        stmp = which(colnames(mm)=="SNR")
        mm = mm[,-stmp]
        mm = cbind(restmp,mm)
        mafout = finalMut[[3]]
        nm = dim(mafout)[1]
        maftmp=apply(resout, 2, function(co) rep(co, each = nm))  ### the matrix is ok, but the data format is not correct
        maftmp = data.frame(maftmp, stringsAsFactors = F)
        maftmp[,2:13] = apply(maftmp[,2:13],2,FUN = function(xx) as.numeric(xx))
        mafout = cbind(maftmp, mafout)
      }else{
        resout = NA
        mafout=NA
      }
      
    }else{
      resout = NA
      mafout=NA
      
    }
    
    return(list(resout,mm,mafout))
  })
  
  flagtmp=sapply(out,FUN=function(outx){
    tmp=outx[[1]][2]
  })
  
  rm(datin)
  gc()
  
  return(out[which(flagtmp>0)])
  
}


chrPeaks=function(chrfile, mafdat, spliting, snrperCut=0.95, ext1side =10000, mutPerKbcut = 5){
  
  dat=read.delim(file = chrfile, header = TRUE, stringsAsFactors = F) 
  
  print(chrfile)
  
  itmp=which(is.infinite(dat$mutrate))
  if(length(itmp)>0){
    dat=dat[-itmp,]
  }
  
  tmp=order(dat$position, dat$mutrate, decreasing = c(FALSE, FALSE))
  odat=dat[tmp,]
  rownames(odat) = as.character(c(1:(dim(odat)[1])))
  
  duppos=which(duplicated(odat$position))
  duppos=unique(odat$position[duppos])
  
  lapply(duppos,FUN=function(apos){
    tmpdat=subset(odat,odat$position==apos)   
    tt=as.numeric(rownames(tmpdat))
    maxo=which.min(tmpdat$mutrate) 
    newdat=tmpdat
    tn=length(tt)
    newmaxo=ceiling(tn/2)
    newdat[newmaxo,]=tmpdat[maxo,]
    if(newmaxo != maxo){
      newdat[maxo,]=tmpdat[newmaxo,]
    } 
    leftdat=NA
    rightdat=NA
    if((newmaxo-1)>=1){
      leftdat=newdat[1:(newmaxo-1),]
      leftdat=leftdat[order(leftdat$mutrate,decreasing = TRUE),]
    }
    if((tn-newmaxo)>=1){
      rightdat=newdat[(newmaxo+1):tn,]
      rightdat=rightdat[order(rightdat$mutrate),]
    }
    outdat=rbind(leftdat,newdat[newmaxo,],rightdat)
    outdat=outdat[!is.na(outdat$mutrate),]
    odat[tt,]=outdat
  })
  
  outfile=strsplit(chrfile, split=spliting)
  outfile=outfile[[1]][2]
  outfile=gsub("\\.","",outfile)
  outfile=gsub("tsv","",outfile)

  datres=peakSearch(datin=odat, snrperc = snrperCut, hKmExt = ext1side, countPerKbcut = mutPerKbcut,
                     chrIn=outfile, mafIn =  mafdat) 
  rm(odat)
  gc()
  
  out1 = do.call(rbind, lapply(datres,FUN=function(xx) xx[[1]]))
  out2 = do.call(rbind, lapply(datres,FUN=function(xx) xx[[2]]))
  out3 = do.call(rbind, lapply(datres,FUN=function(xx) xx[[3]]))
  
  return(list(out1,out2,out3))  
  
}

filterPeaks = function(peakrow, minPeakWidth, maxPeakWidth, noPatCut,mutRateCut){
  tmp1 = ifelse((peakrow[6] - peakrow[5]+1) >= minPeakWidth, 0, 1)
  tmp2 = ifelse((peakrow[2] - peakrow[1] +1) <= maxPeakWidth, 0, 1)
  tmp3 = ifelse(peakrow[3] >= noPatCut, 0, 1)
  tmp4 = ifelse(peakrow[4]>= mutRateCut, 0, 1) 
  tmps = sum(tmp1,tmp2,tmp3,tmp4)
  tmps = sum(tmp1,tmp2,tmp3)
  out = ifelse(tmps>=1, 1, 0)
  return(out)
}

extractMaf1peak=function(awave,mafFileIn){  
  
  names(awave) = c("chromosome", "leftPosition",  "rightPosition")
  awave = data.frame(t(awave), stringsAsFactors = F)
  awave$leftPosition = as.integer(awave$leftPosition)
  awave$rightPosition = as.integer(awave$rightPosition)
  awave = data.table(awave, key =  c("chromosome", "leftPosition",  "rightPosition"))
  
  ##### add one line on July 20, 2018
  ### add one line here on July 20
  mafFileIn = data.table(mafFileIn, key = c("Chromosome", "Start_Position", "End_Position"))
  
  mafsum = foverlaps(x=mafFileIn, y=awave, by.x=c("Chromosome", "Start_Position", "End_Position"), by.y = c("chromosome", "leftPosition",  "rightPosition"))
  mafsum = mafsum[!is.na(mafsum$leftPosition),]
  
  tn=dim(mafsum)[1]
  mutkb=1000 * tn / (awave$rightPosition - awave$leftPosition + 1)
  
  vc=mafsum$Variant_Classification
  vc=table(vc)
  vc=names(which.max(vc))
  
  syms=mafsum$Hugo_Symbo
  syms = syms[which(syms != "Unknown" & syms != "")]
  syms=table(syms)
  syms=names(which.max(syms))
  
  if(length(syms)==0){syms= "Unknown"}
  
  pts = mafsum$Tumor_Sample_Barcode
  pts = length(unique(pts))
  
  outs = list(c(mutkb,pts),c(vc,syms),mafsum)
  return(outs)
  
}


updateWaveouts = function(apeak, mafdat, oldwaves){
  
  finalMut = extractMaf1peak(awave = apeak, mafFileIn = mafdat)
  m1 = as.integer(apeak[2])
  m2 = as.integer(apeak[3])
  mm=subset(oldwaves, oldwaves$LeftPosition >= m1 & oldwaves$RightPosition <= m2)
  mintmp=which.min(mm$peakMin_mutrate)
  ntmp=dim(mm)[1]
  
  meanSNR= mean(mm$SNR)
  
  peakpos=mm$peakMinPosition[mintmp]
  peakmin=mm$peakMin_mutrate[mintmp]
  
  leftpos=m1
  rightpos=m2
  numberP = finalMut[[1]][2]
  
  resout=c(apeak[1],peakpos,peakmin,leftpos,rightpos,numberP,meanSNR,finalMut[[1]][1], finalMut[[2]]) ### 10 items
  resout=data.frame(t(resout), stringsAsFactors = F)
  colnames(resout) = c("chr","peakMinPosition", "peakMin_mutrate", "LeftPosition", "RightPosition","numberOfPatients","meanSNR",
                       "mutPerKb", "mostFreqClassification","mostFreqGeneSymbol")
  resout[,2:8] = apply(resout[,2:8],2,FUN = function(xx) as.numeric(xx))
  mafout = finalMut[[3]]
  nm = dim(mafout)[1]
  maftmp=apply(resout, 2, function(co) rep(co, each = nm))
  maftmp = data.frame(maftmp, stringsAsFactors = F)
  maftmp[,2:8] = apply(maftmp[,2:8],2,FUN = function(xx) as.numeric(xx))
  mafout = cbind(maftmp, mafout)
  
  return(list(resout,mafout))  
  
}


###############################################################
### the following is to apply above functions to get basic wavelet peaks

finalres=lapply(files, chrPeaks, mafdat = maf.full, spliting=sepchr,ext1side = peakExt1side, mutPerKbcut = mutRateCut)
out1 = do.call(rbind, lapply(finalres,FUN=function(xx) xx[[1]]))
out2 = do.call(rbind, lapply(finalres,FUN=function(xx) xx[[2]]))
out3 = do.call(rbind, lapply(finalres,FUN=function(xx) xx[[3]]))

write.table(out1, paste(outname, "waveletMaxChangeSummary.tsv", sep=""),sep = '\t',quote = F) 
write.table(out2, paste(outname, "waveletMaxChangePatientDetail.tsv", sep=""),sep = '\t',quote = F)
write.table(out3, paste(outname, "waveletMaxChangePatientDetailWithMaf.tsv", sep=""),sep = '\t',quote = F)



########################################################################
#### post process steps

## since in the previous step, noPatCut is in %, need to convert to integer based on given number of patients 
allpts = unique(maf.full$Tumor_Sample_Barcode)
allpts = length(allpts)

noPatCut = as.integer( allpts * noPatCut/100)


##################################################################################
##### apply post-process

if(postflag == 1) {
  
  wout1 = out1[,c(4:6,8,9:10)]
  filout1 = apply(wout1, 1, filterPeaks, minPeakWidth, maxPeakWidth, noPatCut,mutRateCut)
  tmp0 = which(filout1<1)
  outf1 = out1[tmp0,]
  
  rangeDat = outf1[,c(1,4:5)]
  dim(rangeDat)  
  
  rangeDat$range = paste(rangeDat$chr, rangeDat$LeftPosition, sep=":")
  rangeDat$range = paste(rangeDat$range, rangeDat$RightPosition, sep="-") 
  
  ##################################################
  ### since some maf and rainstorm chr format might not contain "chr"
  ### should add it in both before using bedr.merge.region function
  tmp = grep("chr", rangeDat$range)
  if(length(tmp) <1){
    rangeDat$range = paste("chr", rangeDat$range, sep="")
  }
  
  ##################################################
  
  mranges = bedr.merge.region(rangeDat$range, distance = distcut) 
  
  mRangeDat = t(sapply(mranges,FUN = function(xx){
    t1=strsplit(xx,split=":")
    t2=strsplit(t1[[1]][2], split="-")
    return(c(t1[[1]][1],t2[[1]][1],t2[[1]][2]))
  }))
  
  ##################################
  ### when tmp add chr, should remove it out again before next step
  if(length(tmp) < 1){
    mRangeDat[,1] = gsub("chr","", mRangeDat[,1])
  }
  
  ##################################
  
  finalout = apply(mRangeDat, 1, updateWaveouts,  mafdat=maf.full, oldwaves = outf1)
  
  sumdat = lapply(finalout, FUN=function(xx) xx[[1]])
  sumdat = do.call(rbind, sumdat)
  
  mafout = lapply(finalout, FUN=function(xx) xx[[2]])
  mafout = do.call(rbind, mafout)
  
  write.table(sumdat, paste(outname, "FilterMergeSummary.tsv", sep=""), sep = '\t',quote = F) 
  write.table(mafout, paste(outname, "FilterMergeDetailWithMaf.tsv", sep=""),sep = '\t',quote = F)
  
}

