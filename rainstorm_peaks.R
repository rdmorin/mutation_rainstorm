#!/usr/bin/env Rscript
#
###########################################
#### R code for peak searching with wavelet
#### Aixiang Jiang, Sep 28, 2017

library(argparse)
library(MassSpecWavelet)
library(maftools)
library(data.table)
#library(parallel)

set.seed(962384) 
#### random seed for peak detection

parser <- ArgumentParser(description="wavelet searching argument");

parser$add_argument(
  "input_files", 
  help="Input path and files", nargs='+', type="character"
);

parser$add_argument(
  "--stringSplit", 
  help="characters before chr# or # in the input file names"
);

parser$add_argument(
  "--input_maf", 
  help="Input maf path and name"
);

parser$add_argument(
  "--output_base_file",  
  help="Output path and name base"
);

parser$add_argument(
           "-patient_minimum", help="minimum number of patients with mutations in a peak for it to be retained", default=4,type="integer");

args = parser$parse_args();

files=args$input_files

sepchr = args$stringSplit
outname = args$output_base_file
maf.file = args$input_maf
patient_minimum = args$patient_minimum
######################################################

peakSearch=function(datin, noPatCut=patient_minimum, snrperc=0.95){  
  newdat=-datin$mutrate
  peakInfo = peakDetectionCWT(newdat) 
  rm(newdat)
  gc()
  
  majorPeakInfo = peakInfo$majorPeakInfo
  res=cbind(majorPeakInfo$peakSNR,majorPeakInfo$allPeakIndex)
  res=res[which(res[,1]>=0),]
  
  #### then, only keep the peaks with 95 percentile of SNR
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
  
  #### calculate median + 25% of mutrate for an entire chr as cutoff
  m25=median(datin$mutrate)+quantile(datin$mutrate, 0.25)
  
  out=apply(ridgeBounds,1,FUN=function(x){
    resout=rep(0, 6)
    tm=NA
    
    m1=max(c((x[1]-12),1))
    m2=min(c((x[2]+12),dim(datin)[1]))
    
    mm=datin[m1:m2,] 
    mintmp=which.min(mm$mutrate)
    ntmp=dim(mm)[1]
 
    tmp1=which(mm$mutrate[1:mintmp]>=m25)
    tmp2=which(mm$mutrate[mintmp:ntmp]>=m25)
    
    rs1=NA
    rs2=NA
   
    if(length(tmp1)>0){
      rs1=max(tmp1)+1
    }else{
      rs1=1
    }
    if(length(tmp2)>0){
      rs2=min(tmp2)-1+mintmp-1
    }else{
      rs2=ntmp
    }
    
    
    if(!is.na(rs1) & !is.na(rs2) & (rs2-rs1+1) >= 4){
      tm=mm[rs1:rs2,]
      tm$SNR=rep(x[3],rs2-rs1+1)
      
      peakpos=mm$position[mintmp]
      peakmin=mm$mutrate[mintmp]
      
      leftpos=min(tm$position)
      rightpos=max(tm$position)
      
      numberP=length(unique(tm$patient))
      
      if(numberP >= noPatCut){
        resout=c(peakpos,peakmin,leftpos,rightpos,numberP,tm$SNR[1])
        pn=dim(tm)[1]
        tm=cbind(rep(peakpos,pn),rep(leftpos,pn),rep(rightpos,pn),tm)
      }
      #}
    }
    
    return(list(resout,tm))
  })
  
  flagtmp=sapply(out,FUN=function(outx){
    tmp=outx[[1]][1]
    return(sum(tmp))
  })
  
  rm(datin)
  gc()
  
  return(out[which(flagtmp>0)])
  
}

##### end of the new peak searching function 
###############################################


#############################################
chrPeaks=function(chrfile, spliting, snrperCut=0.95){
  
  dat=read.table(file = chrfile, sep = '\t', header = TRUE, stringsAsFactors = F) 
  
  #### to avoid problems, remove -inf and NA rows
  itmp=which(is.infinite(dat$mutrate))
  if(length(itmp)>0){
    dat=dat[-itmp,]
  }
  dat=dat[!is.na(dat$mutrate),]
  
  tmp=order(dat$position, dat$mutrate, decreasing = c(FALSE, FALSE))
  odat=dat[tmp,]
  rownames(odat) = as.character(c(1:(dim(odat)[1])))
  
  duppos=which(duplicated(odat$position))
  duppos=unique(odat$position[duppos])

  lapply(duppos,FUN=function(apos){
    tmpdat=subset(odat,odat$position==apos)   #### treat odat as a global obj without passing it into this function
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
  
  datres=peakSearch(datin=odat, snrperc = snrperCut)
  
  rm(odat)
  gc()
  
  numres=t(sapply(datres,FUN=function(x){
    return(x[[1]])
  }))
  colnames(numres)=c("peakPosition","peakMinValue","leftPosition","rightPosition","numberPatients","SNR")

  numres=data.frame(numres)
  outfile=strsplit(chrfile, split=spliting)
  outfile=outfile[[1]][2]
  outfile=gsub("\\.","",outfile)
  outfile=gsub("tsv","",outfile)
  numres$chromosome=rep(outfile,dim(numres)[1])
  
  patres=lapply(datres,FUN=function(x){
    x=x[[2]]
    x$chromosome=rep(outfile,dim(x)[1])
    colnames(x)[1:3]=c("peakPosition","leftPosition","rightPosition")
    return(x)
  })
  
  allranges=numres[,3:4]
  allouts=apply(allranges,1,FUN=function(x){
    tmp=apply(allranges,1,FUN=function(oneRange){
      if(oneRange[1]>x[2] | oneRange[2]<x[1]){
        return(0)
      }else{
        return(1)
      }
    })
    tmp=which(tmp>0)
    
    if(length(tmp)>1){
      tmpdat=numres[tmp,]
  
      numtmp=tmpdat[1,]
      numtmp[1]=tmpdat[which.min(tmpdat[,2]),1]
      numtmp[2]=min(tmpdat[,2])
      numtmp[3]=min(tmpdat[,3])
      numtmp[4]=max(tmpdat[,4])
      numtmp[6]=min(tmpdat[,6])
      pattmp=do.call(rbind, lapply(tmp,FUN=function(yy){patres[[yy]]}))
      pattmp=pattmp[!duplicated(pattmp[,c(6:7)]),]
      ntmp=dim(pattmp)[1]
      numtmp[5]=length(unique(pattmp$patient))
      pattmp$peakPosition=rep(as.numeric(numtmp[1]),ntmp)
      pattmp$leftPosition=rep(as.numeric(numtmp[3]),ntmp)
      pattmp$rightPosition=rep(as.numeric(numtmp[4]),ntmp)
      pattmp$SNR = rep(as.numeric(numtmp[6]), ntmp)
    }else{
      numtmp=numres[tmp,]
      pattmp=patres[[tmp]]
    }
    
    ##### add extra columns for mean mutrate, sd mutrate, mut_per_kb, 
    ntmp=dim(pattmp)[1]
    numtmp[8]=mean(pattmp$mutrate)
    numtmp[9]=sd(pattmp$mutrate)
    numtmp[10]=1000 * ntmp / (numtmp[4] - numtmp[3]  + 1)
    pattmp$meanOfmutrate=rep(as.numeric(numtmp[8]),ntmp)
    pattmp$sdOfmutrate=rep(as.numeric(numtmp[9]),ntmp)
    pattmp$mutPerKb=rep(as.numeric(numtmp[10]),ntmp)
    
    return(list(numtmp, pattmp))
    
  })
  
  rm(datres)
  gc()
  
  
  #### build up two data.frame to return
  allnum=sapply(allouts,FUN=function(reslist){
    return(reslist[1])
  })
  allnum=do.call(rbind,allnum)
  colnames(allnum)[8:10]=c("meanOfmutrate","sdOfmutrate","mutPerKb")
  
  allpat=sapply(allouts,FUN=function(reslist){
    return(reslist[2])
  })
  allpat=do.call(rbind,allpat)
  
  ##### remove duplicated
  allnum=allnum[!duplicated(allnum$peakPosition),]
  allpat=allpat[!duplicated(allpat[,c("patient","position")]),]
  
  rm(allouts)
  gc()
  
  return(list(allnum,allpat))
  
}

finalres=sapply(files, chrPeaks, spliting=sepchr)

##### build up two matrix across chrs
nn=length(finalres)
ii=c(1:nn)
m1=which(ii %% 2 == 1)
m2=which(ii %% 2 == 0)

m1res=finalres[m1]
m1res=do.call(rbind,m1res)
m2res=finalres[m2]
m2res=do.call(rbind,m2res)

write.table(m1res, paste(outname, "waveletSummary.tsv", sep=""),sep = '\t',quote = F, row.names = F)
write.table(m2res, paste(outname, "waveletPatientDetail.tsv", sep=""),sep = '\t',quote = F, row.names = F)

vc = c("3'Flank","3'UTR","5'Flank","5'UTR","Frame_Shift_Del","Frame_Shift_Ins","IGR","In_Frame_Del","In_Frame_Ins","Intron","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","RNA","Silent","Splice_Region","Splice_Site","Translation_Start_Site")

maf.full = read.maf(maf.file,useAll = T, vc_nonSyn=vc)
mafs=data.frame(maf.full@data)
rm(maf.full)
gc()

withMaf=function(awave,patRains, mafFile, mutrateKbCut){
  #### calculate mutrate per kb, if it pass a cutoff, then move on, otherwise return NA
  tmp=subset(mafFile,mafFile$Chromosome == as.character(awave[7]) & mafFile$Start_Position >= as.integer(awave[3]) 
             & mafFile$Start_Position <= as.integer(awave[4]))
  
  tn=dim(tmp)[1]
  mutkb=1000 * tn / ( as.integer(awave[4]) -  as.integer(awave[3])  + 1)
  outs=NA
  if(mutkb>=mutrateKbCut){
    #### find the most freq variant_classification, since each entry in mafFile has an efficient Variant_Classification
    #### we do not need to worry about NA or null...
    vc=tmp$Variant_Classification
    vc=table(vc)
    vc=names(which.max(vc))
    
    #### merge thing, no need the chromosome within a peak range, but do need position and patient
    #### find the subset from the patRains 
    ptmp=subset(patRains, patRains$chromosome==as.character(awave[7]) & patRains$peakPosition == as.integer(awave[1]))
    #### for tmp, we need Start_Position and Tumor_Sample_Barcode
    #### for ptmp, need position and patient
    alltmp=merge(ptmp, tmp, by.x=c("position","patient"), by.y=c("Start_Position","Tumor_Sample_Barcode"), all=TRUE)
    
    ####for these have NA, input: peakPosition leftPosition rightPosition SNR chromosome meanOfmutrate sdOfmutrate mutPerKb
    tt=which(is.na(alltmp$peakPosition))
    ttt=alltmp[-tt,]
    alltmp[tt,c(3:5,9:13)]=ttt[1,c(3:5,9:13)]
    
    nall=dim(alltmp)[1]
    #### add more columns in, awave[2], new number of patient column, mutkb, and vc
    
    newcol=data.frame(cbind(rep(as.numeric(awave[2]),nall),rep(length(unique(alltmp$patient)),nall),rep(mutkb, nall)))
    colnames(newcol)=c("peakMinValue","numberPatients","mutPerKbMaf")
    newcol$mostFreqVClassification = rep(vc,nall)
    
    outs=data.frame(cbind(alltmp[,c(10,4:5,3)], newcol,alltmp[,-c(10,4:5,3)]))
  }
  
  return(outs)
}

mafm1res=apply(m1res,1,withMaf,patRains=m2res, mafFile=mafs, mutrateKbCut=6)
mafm1res=do.call(rbind, mafm1res)
#### remove all of the NAs
mafm1res=mafm1res[!is.na(mafm1res$chromosome),]

write.table(mafm1res,paste(outname,"waveletPatientDetail_withMaf.tsv", sep=""),sep = '\t',quote = F, row.names = F)

#### get the summary one row per peak range as well
sumres=mafm1res[,c(1:8,14:19)]

tmp=unique(sumres[,c("chromosome", "leftPosition", "rightPosition", "peakPosition")])
getHGall=function(x){
  subdat=subset(sumres, sumres$chromosome==as.character(x[1]) & sumres$leftPosition == as.integer(x[2])
                & sumres$rightPosition == as.integer(x[3]) & sumres$peakPosition == as.integer(x[4]))
  hg=table(subdat$Hugo_Symbol)
  tt=NA
  if(length(hg)>0){
    tt=subdat$Hugo_Symbol[which.max(hg)]
  }
  subdat$Hugo_Symbol[1]=tt
  subdat=subdat[1,]
}
ss=apply(tmp,1,getHGall)
ss=do.call(rbind,ss)


write.table(ss,paste(outname, "waveletSummary_withMaf.tsv",sep=""),sep = '\t',quote = F, row.names = F)


gc()

quit()



