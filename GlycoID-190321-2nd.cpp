#include <Rcpp.h> 
#include   <stdlib.h>
#include <math.h>  
using namespace Rcpp; 
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <map>
#include <algorithm>
#include <iterator>
#include<vector>

// [[Rcpp::export]]

NumericVector deisotope(NumericVector x, NumericVector y,NumericVector z){
  
  
  int m=x.size();
  NumericVector iso(m);
  std::fill(iso.begin(), iso.end(), 0);
  
  
  for(int i=0;i<m;i++)
  {
    int minsize=std::min(i+5,m-1);
    for(int j=i+1;j<minsize;j++)
    {
      if((y[i]>350)&&((fabs(y[j]-y[i]-1.0034) <0.005) |(fabs(y[j]-y[i]-0.50017)<0.0025))&&(x[i]==x[j])&&(z[i]>2000))
        iso[j]=1;  
    }
  }
  
  for(int i=0;i<m;i++)
  {
    int minsize=std::min(i+5,m-1);
    for(int j=i+1;j<minsize;j++)
    {
      if((y[i]>350)&&((fabs(y[j]-y[i]-1.0034) <0.005) |(fabs(y[j]-y[i]-0.50017)<0.0025))&&(x[i]==x[j])&&(z[i]>2000))
        iso[j]=1;  
    }
  }
  
  return wrap(iso);
  
}

// [[Rcpp::export]]

NumericVector  calmass(std::vector< std::string > seq1){
  
  
  int m=seq1.size();
  
  std::vector< double > pepmass;
  
  
  std::map<char,double> AAmass;
  AAmass.insert ( std::pair<char,double>('A',71.0371) );
  AAmass.insert ( std::pair<char,double>('R',156.1011) );
  AAmass.insert ( std::pair<char,double>('N',114.0429) );
  AAmass.insert ( std::pair<char,double>('D',115.0269) );
  AAmass.insert ( std::pair<char,double>('C',160.0338) );
  AAmass.insert ( std::pair<char,double>('E',129.0426) );
  AAmass.insert ( std::pair<char,double>('Q',128.0586) );
  AAmass.insert ( std::pair<char,double>('G',57.0215) );
  AAmass.insert ( std::pair<char,double>('H',137.0589) );
  AAmass.insert ( std::pair<char,double>('I',113.0841) );
  AAmass.insert ( std::pair<char,double>('L',113.0841) );
  AAmass.insert ( std::pair<char,double>('K',128.0950) );
  AAmass.insert ( std::pair<char,double>('M',131.0405) );
  AAmass.insert ( std::pair<char,double>('F',147.0684) );
  AAmass.insert ( std::pair<char,double>('P',97.0528) );
  AAmass.insert ( std::pair<char,double>('S',87.0320) );
  AAmass.insert ( std::pair<char,double>('T',101.0477) );
  AAmass.insert ( std::pair<char,double>('W',186.0793) );
  AAmass.insert ( std::pair<char,double>('Y',163.0633) );
  AAmass.insert ( std::pair<char,double>('V',99.0684) );
  
  
  for(int i=0;i<m;i++)
  {
    std::string pepsequen1=seq1[i];
    double sumass=0;
    
    for (std::string::iterator it=pepsequen1.begin(); it!=pepsequen1.end(); ++it)
    {
      double residulmass= AAmass.find(*it)->second;
      sumass=sumass+residulmass;
    }
    
    sumass=sumass+1.0078+18.0106;
    pepmass.push_back(sumass);
    
    
  }
  
  return wrap(pepmass);
}


// [[Rcpp::export]]

int mass300pos(NumericVector x){
  int m=x.size();
  int i;
  for(i=0;i<m;i++){if(x[i]>250) break;}
  
  
  
  return(i);
  
}

// [[Rcpp::export]]

int spectragap(NumericVector x,NumericVector y){
  
  int m=x.size();
  int z=0;
  
  
  for(int i=m;i>1;i--)
  {
    int ion1=x[i];
    int msize=std::max(i-11,0);
    for(int j=i-1;j>msize-1;j--)
    {
      int tempg=ion1-x[j];
      bool gaptrue=is_true(any(abs(tempg-y)<0.06));
      if(gaptrue)   z++;
    }
    
  }
  
  return z;
  
}

// [[Rcpp::export]]

DataFrame getspectraintensity(NumericVector x, NumericVector y,NumericVector z,int t){
  
  NumericVector mz;
  NumericVector intensity;
  int m=x.size();
  
  
  for(int i=0;i<m;i++)
  {
    if ( x[i]<t) 
    {
      continue;
    }
    
    else if ( x[i]==t) 
    {
      mz.push_back(y[i]);
      intensity.push_back(z[i]);
      
    }
    else
      break;
  }
  
  return Rcpp::DataFrame::create(Named("mz")=mz,Named("intensity")=intensity);
  
}
// [[Rcpp::export]]

NumericVector getspectra(NumericVector x, NumericVector y,int t){
  
  NumericVector z;
  int m=x.size();
  
  
 for(int i=0;i<m;i++)
 {
   if ( x[i]<t) 
   {
     continue;
   }
   
   else if ( x[i]==t) 
   {
     z.push_back(y[i]);
   }
   else
     break;
 }
     
  return wrap(z);
  
}


// [[Rcpp::export]]

int spectracomp(NumericVector x, NumericVector y,double delta){
  
  int m = x.size();
  int n = y.size(); 
  int t=0;
  int i=0;
  int j=0;
  
  while(i<m && j<n)
  {
    
    if ((x[i]-y[j])<delta&&(x[i]-y[j])>(-delta))
    {
      i=i+1;
      j=j+1;
      t++;
    }
    else if(x[i]<=y[j]-delta) i=i+1;
    else if(x[i]>=y[j]+delta) j=j+1;
  }
  
  return t;
}

// [[Rcpp::export]]

DataFrame Tspectracomp(NumericVector mgfindex1,NumericVector pepmz1,NumericVector Rtime1,NumericVector charge1,
                       NumericVector mgfindex2,NumericVector pepmz2,NumericVector Rtime2,NumericVector charge2,
                       NumericVector mgfindex11,NumericVector mz1,NumericVector mgfindex12,NumericVector mz2) {
  
  
  int start=0;
  int m = mgfindex1.size();
  int n = mgfindex2.size();
  std::vector<int> score;
  std::vector<int> scoreoxonium;
  std::vector<int> gapionnum;
  std::vector<int> mgfi;
  std::vector<int> mgfj;
  std::vector<int> sizei;
  std::vector<int> sizej;
  std::vector<double> deltamass;
  
  NumericVector spectraoxonmium = NumericVector::create(
    126.0550,138.0550,145.0500,163.0606,168.066,186.0166,
    204.0872,274.0926,292.1032,325.1134,366.14);
  
  NumericVector gapions = NumericVector::create(
   66.0212,73.0259,81.0264,88.0161, 101.5397,132.0423,145.5477,
   146.0579,153.5451,162.0528,176.0321,203.0794,291.0954,307.0903);
  
  for(int i=0;i<m-1;i++)
  {
    NumericVector spectra1=getspectra(mgfindex11,mz1,mgfindex1[i]);
    int scoretempoxo=spectracomp(spectra1,spectraoxonmium,0.02);
    int tempgapions=spectragap(spectra1,gapions);
    double pepmassi=pepmz1[i]*charge1[i];
    
    int pos3001=mass300pos(spectra1);//**mass below 300 are filtered
    
    spectra1.erase(spectra1.begin(),spectra1.begin()+pos3001);//**mass below 300 are filtered
    
    
    
    for(int j=start;j<n;j++)
    {
      if((Rtime2[j]-Rtime1[i])<=(-240))
      {
        start++;
        
      }
      else if((-240<(Rtime2[j]-Rtime1[i]))&&((Rtime2[j]-Rtime1[i])<240))
      {
        double pepmassj=pepmz2[j]*charge2[j];
        double tempdeltamass=pepmassi-pepmassj;
      
        NumericVector spectra2=getspectra(mgfindex12,mz2,mgfindex2[j]);
        int pos3002=mass300pos(spectra2);//**mass below 300 are filtered
        spectra2.erase(spectra2.begin(),spectra2.begin()+pos3002);//**mass below 300 are filtered
        int scoretemp=spectracomp(spectra1,spectra2,0.05);
        if(scoretemp>5){
        score.push_back(scoretemp);
        scoreoxonium.push_back(scoretempoxo);
        gapionnum.push_back(tempgapions);
        deltamass.push_back(tempdeltamass);
        mgfi.push_back(mgfindex1[i]);
        mgfj.push_back(mgfindex2[j]);
        sizei.push_back(spectra1.size());
        sizej.push_back(spectra2.size());
        }
      }
      else continue;
  }
  }
  return Rcpp::DataFrame::create(Named("mgfindex1")=mgfi,Named("mgfindex2")=mgfj,Named("size1")=sizei,Named("size2")=sizej,
                                       Named("score")=score,Named("deltamass")=deltamass,Named("scoreoxonium")=scoreoxonium, Named("gapionsnum")=gapionnum);
}


// [[Rcpp::export]]

DataFrame search3rd(NumericVector mgfindex1,NumericVector scanindex1,NumericVector scanindex2,NumericVector probscore1,NumericVector simliarity1,
                    std::vector< std::string > peptideseq,NumericVector mgfindex11,NumericVector mz1,
                    std::vector< std::string > sequence,NumericVector peplistmass) {
  
  
  
  int m = mgfindex1.size();
  int n = peplistmass.size();
  std::vector<int> scani1;
  std::vector<int> scani2;
  std::vector<double> probscore;
  std::vector<double> simliarity;
  std::vector<int> mgfi; 
  std::vector< std::string > pepseq;
  std::vector<int> Nprecur;
  
  for(int i=0;i<m-1;i++)
  {
    NumericVector spectra=getspectra(mgfindex11,mz1,mgfindex1[i]+1);
    
    std::string pepi=peptideseq[i];
    
    for(int j=0;j<n;j++)
    {
      std::string pepj=sequence[j];
      int sfind1=pepi.find(pepj);
      int sfind2=pepj.find(pepi);
      
      if((sfind1!=(-1))|(sfind2!=(-1))){
        
        
        NumericVector precursorlist=NumericVector::create((peplistmass[j]+1.0078)/2,(peplistmass[j]+1.0078+203.0794)/2,peplistmass[j],peplistmass[j]+203.0794);
        int massmatch=spectracomp(spectra,precursorlist,0.05);
        
        mgfi.push_back(mgfindex1[i]); 
        scani1.push_back(scanindex1[i]);
        scani2.push_back(scanindex2[i]);
        probscore.push_back(probscore1[i]);
        simliarity.push_back(simliarity1[i]);
        pepseq.push_back(pepj);
        Nprecur.push_back(massmatch);
        
      }
      
    }
  }
  return Rcpp::DataFrame::create(Named("mgfindex1")=mgfi,Named("scanindex1")=scani1,Named("scanindex2")=scani2,Named("probscore")=probscore,Named("simliarity")=simliarity,
                                 Named("pepseq")=pepseq,Named("Nprecur")=Nprecur);
}

/*** R

setwd("D:\\Rfiles\\Glycan")
library(data.table)
library(Rcpp)
library(inline)
library(ggplot2)
  

spectra<-function(df)
{
  sax<-df[,1:2]
  colnames(sax)<-c("mz","intensity")
  pepmassindex<-which(sax$mz=="TITLE")
  pepmass<-sax$intensity[pepmassindex+1]
  pepmass<-strsplit(pepmass," ")
  pepmass<-sapply(pepmass,function(x) x[1])
  charge<-sax$intensity[pepmassindex+2]
  charge<-strsplit(charge,NULL)
  charge<-sapply(charge,function(x) x[1])
  Rtime<-sax$intensity[pepmassindex+3]
  scanindex<-sax$intensity[pepmassindex+4]
  mgfindex=seq(1:length(pepmass))
  parention<-data.frame(mgfindex,pepmass,charge,Rtime,scanindex)
  remove_index<-c(pepmassindex,pepmassindex+1,pepmassindex+2,pepmassindex+3,pepmassindex+4,pepmassindex+5)
  fragmention<-sax$mz[-remove_index]
  index<-which((fragmention=="BEGIN IONS")|(fragmention=="\nBEGIN IONS"))
  cumdex<-rep(0,times=length(fragmention))
  cumdex[index]<-1
  cumdex[1]<-1
  mgfindex<-cumsum(cumdex)
  fragmention<-as.character(fragmention)
  fragmention<-strsplit(fragmention," ")
  mz<-sapply(fragmention,function(x) x[1])
  intensity<-sapply(fragmention,function(x) x[2])
  fragmention<-data.frame(mgfindex,mz,intensity)
  index<-which(fragmention$intensity=="IONS")
  fragmention<-fragmention[-index,]
  rm(sax)
  return (list("parention"=parention,"fragmention"=fragmention))
}











sax<-fread("guoyu-C+F-3.mgf",sep="=", fill=TRUE,blank.lines.skip=TRUE,data.table=TRUE)
  sax<-data.frame(sax)
  
  
  sax<-spectra(sax)
  parention1<-sax$parention
  spectralist1<-sax$fragmention
  write.csv(spectralist1,file="spectralist1.csv")
  rm(spectralist1)
  spectralist1<-fread("spectralist1.csv")
  
  
  
  
  
  spectralist1<-as.data.table(spectralist1)
  spectralist1<-spectralist1[,.SD[1],by=.(mgfindex,mz)]
setkey(spectralist1,mgfindex,intensity)
  spectralist1<-spectralist1[intensity>2000,]
spectralist1<-spectralist1[,tail(.SD,150),by=mgfindex]
mgfindex<-as.numeric(spectralist1$mgfindex)
  mz<-as.numeric(spectralist1$mz)
  spectralist1<-spectralist1[order(mgfindex,mz),]
write.csv(parention1,file="parention1.csv")
  rm(parention1)
  parention1<-read.csv("parention1.csv")
  mgfindex1<-as.numeric(parention1$mgfindex)
  pepmz1<-parention1$pepmass
  Rtime1<-parention1$Rtime
  charge1<-as.numeric(parention1$charge)
  
  
  
  spectralist1mgf<-as.numeric(spectralist1$mgfindex)
  spectralist1mz<-as.numeric(spectralist1$mz)
  spectralist1intensity<-as.numeric(spectralist1$intensity)
  isotope<-deisotope(spectralist1mgf,spectralist1mz,spectralist1intensity)
  spectralisttemp<-data.frame(spectralist1,isotope)
  index<-which(spectralisttemp$isotope==0)
  spectralist1<-spectralisttemp[index,]
  spectralist1<-spectralist1[,1:4]
  mgfindex11<-as.numeric(spectralist1$mgfindex)
  mz1<-as.numeric(spectralist1$mz)
  intensity1<-as.numeric(spectralist1$intensity)
  
  

    
    
  psms<-read.csv("guoyu-C+F-ID-psms.csv")
  psms<-as.data.table(psms)
  
  psms<-psms[,.(pepsequen,scanindex)]

  spectraref<-merge(parention1,psms,by="scanindex")
  spectraref<-na.omit(spectraref)
  spectraref<-as.data.table(spectraref)
  mgfindex<-spectraref$mgfindex
  spectraref<-spectraref[order(mgfindex),]
  spectraref<-as.data.table(spectraref)
  spectraref2<-spectraref[,.(mgfindex,scanindex)]
  write.csv(spectraref,file="spectraref.csv")
  rm(spectraref)
  spectraref<-read.csv("spectraref.csv")


  spectralist2<-spectralist1
  spectralist2<-merge(spectralist2,spectraref2,by="mgfindex")
  spectralist2<-na.omit(spectralist2)
 
  mgfindex<-as.numeric(spectralist2$mgfindex)
  mz<-as.numeric(spectralist2$mz)
  spectralist2<-spectralist2[order(mgfindex,mz),]
  mgfindex2<-as.numeric(spectraref$mgfindex)
  pepmz2<-spectraref$pepmass
  Rtime2<-spectraref$Rtime
  charge2<-as.numeric(spectraref$charge)
  mgfindex22<-as.numeric(spectralist2$mgfindex)
  mz2<-as.numeric(spectralist2$mz)
  intensity2<-as.numeric(spectralist2$intensity)
  rm(sax)
  rm(spectralist2)
  comparedspectra1<-data.frame(mgfindex1,pepmz1,Rtime1,charge1)
  comparedspectra2<-data.frame(mgfindex2,pepmz2,Rtime2,charge2)


  



  
score<-Tspectracomp(mgfindex1,pepmz1,Rtime1,charge1,mgfindex2,pepmz2,Rtime2,charge2,
                      mgfindex11,mz1,mgfindex22,mz2)
index<-which(score$score>12)
score<-score[index,]
index<-which(score$scoreoxonium>1)
score<-score[index,]


leng<-length(score$mgfindex1)
scanindex1<-rep(NA,times=leng)
scanindex2<-rep(NA,times=leng)
peptideseq<-rep(NA,times=leng)
pre_mz1<-rep(NA,times=leng)
pre_charge1<-rep(NA,times=leng)
pre_mz2<-rep(NA,times=leng)
pre_charge2<-rep(NA,times=leng)
Rtimec1<-rep(NA,times=leng)
Rtimec2<-rep(NA,times=leng)
for(i in 1:leng)
{
  mgfindex1temp<-score$mgfindex1[i]
  indextemp1<-which(parention1$mgfindex==mgfindex1temp)
  scanindex1[i]<-parention1$scanindex[indextemp1]
  pre_mz1[i]<-parention1$pepmass[indextemp1]
  pre_charge1[i]<-parention1$charge[indextemp1]
  Rtimec1[i]<-parention1$Rtime[indextemp1]
  mgfindex2temp<-score$mgfindex2[i]
  indextemp2<-which(parention1$mgfindex==mgfindex2temp)
  scanindex2[i]<-parention1$scanindex[indextemp2]
  
  indextemp2<-which(spectraref$mgfindex==mgfindex2temp)
  peptideseq[i]<-as.character(spectraref$pepsequen[indextemp2])
  pre_mz2[i]<-spectraref$pepmass[indextemp2]
  pre_charge2[i]<-spectraref$charge[indextemp2]
  Rtimec2[i]<-spectraref$Rtime[indextemp2]
}

score<-data.frame(score,scanindex1,scanindex2,peptideseq,pre_mz1,pre_charge1,pre_mz2,pre_charge2,Rtimec1,Rtimec2)

probscore<-1-ppois(score$score,(score$size1*score$size2/1000))
probscore<-(-7.2)*log10(probscore)


newspectra<-data.frame(score,probscore)



leng<-length(newspectra$mgfindex1)
simlarity<-rep(NA,times=leng)
commonions<-rep(NA,times=leng)
spectraoxo<-c(126.0550,138.0550,145.0500,163.0606,168.066,186.0166,
              204.0872,274.0926,292.1032,325.1134,366.14)

for(i in 1:leng)
{
  mgfindex1temp<-newspectra$mgfindex1[i]
  spectra1<-getspectraintensity(mgfindex11,mz1,intensity1,mgfindex1temp)
  colnames(spectra1)<-c("mz1","intensity1")
  spectra1<-as.data.table(spectra1)
  spectra1<-spectra1[,start:=mz1-0.03]
  spectra1<-spectra1[,end:=mz1+0.03]
  setkey(spectra1,start,end)
  
  mgfindex2temp<-newspectra$mgfindex2[i]
  spectra2<-getspectraintensity(mgfindex22,mz2,intensity2,mgfindex2temp)

  colnames(spectra2)<-c("mz2","intensity2")
  spectra2<-as.data.table(spectra2)
  spectra2<-spectra2[,start:=mz2-0.03]
  spectra2<-spectra2[,end:=mz2+0.03]
  setkey(spectra2,start,end)
  commonspectra<-foverlaps(spectra1,spectra2,mult="first")
  index<-which(!is.na(commonspectra$mz2))
  commonspectra<-commonspectra[index,]
  commonions[i]<-length(commonspectra$mz2)
  
  intensity11<-commonspectra$intensity1
  intensity22<-commonspectra$intensity2
  simlarity[i]<-sum(intensity11*intensity22)/sqrt((sum(intensity11^2))*(sum(intensity22^2)))
 }

newspectra<-data.frame(newspectra,commonions,simlarity)

leng<-length(newspectra$mgfindex1)
oxoionsDT<-rep(NA,times=leng)
spectraoxo<-c(126.0550,138.0550,145.0500,163.0606,168.066,186.0166,
                204.0872,274.0926,292.1032,325.1134,366.14)
spectraoxo<-data.frame(spectraoxo)
colnames(spectraoxo)<-c("mzoxo")
  
  for(i in 1:leng)
  {
    mgfindex1temp<-newspectra$mgfindex1[i]
    spectra1<-getspectraintensity(mgfindex11,mz1,intensity1,mgfindex1temp)
    colnames(spectra1)<-c("mz1","intensity1")
    spectra1<-as.data.table(spectra1)
    spectra1<-spectra1[,start:=mz1-0.02]
    spectra1<-spectra1[,end:=mz1+0.02]
    setkey(spectra1,start,end)
    
    spectra1trun<-spectra1[,intenratio:=intensity1/max(intensity1)]
    spectra1trun<-spectra1trun[intenratio>0.25]
    spectra1trun<-spectra1trun[mz1<400,]
    setkey(spectra1trun,start,end)
    
    spectraoxo<-as.data.table(spectraoxo)
    spectraoxo<-spectraoxo[,start:=mzoxo-0.001]
    spectraoxo<-spectraoxo[,end:=mzoxo+0.001]
    setkey(spectraoxo,start,end)
    
    oxocommomspectra<-foverlaps(spectraoxo,spectra1trun,mult="first")
    index<-which(!is.na(oxocommomspectra$mz1))
    oxocommomspectra<-oxocommomspectra[index,]
    oxoionsDT[i]<-length(oxocommomspectra$mz1)
  }
pepmass<-calmass(as.character(newspectra$peptideseq))
newspectra<-data.frame(newspectra,pepmass)
deltamss<-newspectra$pre_mz1*newspectra$pre_charge1-1.0078*(newspectra$pre_charge1-1)-pepmass
newspectra<-data.frame(newspectra,deltamss)

index<-which(newspectra$probscore>20)

mgfindex1<-newspectra$mgfindex1
peptideseq<-as.character(newspectra$peptideseq)
    scanindex1<-newspectra$scanindex1
    scanindex2<-newspectra$scanindex2
    index<-which(newspectra$probscore>200)
    newspectra$probscore[index]=200
  probscore<-newspectra$probscore
    simlarity<-newspectra$simlarity
    
    peptidelist<-read.csv("guoyu-PM1-list.csv")
    peptidelist<-as.character(peptidelist$Sequence)
    pepmass<-calmass(peptidelist)
    newresult<-search3rd(mgfindex1,scanindex1,scanindex2,probscore,simlarity,peptideseq,mgfindex11,mz1,peptidelist,pepmass)
    
    index<-which(newresult$Nprecur>0)
    score<-newresult[index,]
  
  
  leng<-length(score$mgfindex1)
    
    pre_mz1<-rep(NA,times=leng)
    pre_charge1<-rep(NA,times=leng)
    Rtimec1<-rep(NA,times=leng)
    Rtimec2<-rep(NA,times=leng)
    for(i in 1:leng)
    {
      mgfindex1temp<-score$mgfindex1[i]
      indextemp1<-which(newspectra$mgfindex1==mgfindex1temp)
      pre_mz1[i]<-newspectra$pre_mz1[indextemp1[1]]
      pre_charge1[i]<-newspectra$pre_charge1[indextemp1[1]]
      Rtimec1[i]<-newspectra$Rtimec1[indextemp1[1]]
    }
    
    
    newspectra2<-data.frame(score,pre_mz1,pre_charge1,Rtimec1)
      
      write.csv(newspectra,file="score.csv")
*/
