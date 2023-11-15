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

NumericVector  singlespectragen(std::string pepsequen1){
  
  std::vector< double > mz;
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
  
  
  double sumass=1.0078;
  for (std::string::iterator it=pepsequen1.begin(); it!=pepsequen1.end()-1; ++it)
  {
    double residulmass= AAmass.find(*it)->second;
    sumass=sumass+residulmass;
    mz.push_back(sumass);
  }
  
  
  sumass=19.0184;
  for (std::string::iterator it=pepsequen1.end()-1; it!=pepsequen1.begin(); --it)
  {
    double residulmass= AAmass.find(*it)->second;
    sumass=sumass+residulmass;
    mz.push_back(sumass);
  }
  
  std::sort(mz.begin(),mz.end());
  return wrap(mz);
}



// [[Rcpp::export]]

NumericVector vectorsort(NumericVector x){
  NumericVector spectra;
  
  for (int i=0; i<x.size(); i++){
    double t=x[i];
    if(t<1800&&t>71){spectra.push_back(t);}
    else if(t>1600&&t<3600){
      t=(t+1.0078)/2;
      spectra.push_back(t);}
    }
  
  NumericVector spectra2;
  spectra2=unique(spectra);
  std::sort(spectra2.begin(), spectra2.end());
  
  return(spectra2);
  
}

// [[Rcpp::export]]

NumericVector vectorsortcz(NumericVector x, NumericVector y){
  NumericVector spectra;
  
  for (int i=0; i<x.size(); i++){
    double t=x[i];
    if(t>1600&&t<3600){
      t=(t+1.0078)/2;
      spectra.push_back(t);}
    else if(t<1600&&t>1000)
      {
      spectra.push_back(t);
      t=(t+1.0078)/2;
      spectra.push_back(t);
      }
    else if(t<1000){spectra.push_back(t);}
   }
  
  for (int i=0; i<y.size(); i++){
    double t=y[i];
    if(t>1600&&t<3600){
      t=(t+1.0078)/2;
      spectra.push_back(t);}
    else if(t<1600&&t>1000)
    {
      spectra.push_back(t);
      t=(t+1.0078)/2;
      spectra.push_back(t);
    }
    else if(t<1000){spectra.push_back(t);}
  }
    
  NumericVector spectra2=unique(spectra);
  std::sort(spectra2.begin(), spectra2.end());
  
  
  
  return(spectra2);
  
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
  for(i=0;i<m;i++){if(x[i]>200) break;}
  
  
  
  return(i);
  
}

// [[Rcpp::export]]

NumericVector calz(NumericVector x, NumericVector y,NumericVector z){
  
  
  int m=x.size();
  NumericVector iso(m);
  std::fill(iso.begin(), iso.end(), 1);
  
  
  for(int i=0;i<m;i++)
  {
    int minsize=std::min(i+5,m-1);
    for(int j=i+1;j<minsize;j++)
    {
      if((y[i]>250)&&((fabs(y[j]-y[i]-0.50017)<0.0025))&&(x[i]==x[j])&&(z[i]>2000))
        iso[i]=2;  
    }
  }
  
  for(int i=0;i<m;i++)
  {
    int minsize=std::min(i+5,m-1);
    for(int j=i+1;j<minsize;j++)
    {
      if((y[i]>250)&&((fabs(y[j]-y[i]-0.3344) <0.005) )&&(x[i]==x[j])&&(z[i]>2000))
        iso[i]=3;  
    }
  }
  
  return wrap(iso);
  
}



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

int spectragap(NumericVector x,NumericVector y){
  
  int m=x.size();
  int z=0;
  
  
  for(int i=m;i>1;i--)
  {
    double ion1=x[i];
    int msize=std::max(i-15,0);
    for(int j=i-1;j>msize-1;j--)
    {
      double tempg=ion1-x[j];
      bool gaptrue=is_true(any(abs(tempg-y)<0.02));
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

DataFrame Tspectracomp(NumericVector mgfindex1,NumericVector pepmz1,IntegerVector Rtime1,NumericVector charge1,
                        NumericVector mgfindexexp,NumericVector mzexp,NumericVector pepmass2,IntegerVector scanindex,NumericVector Rtime2,std::vector< std::string > sequence) {
  
  
 
  int m = mgfindex1.size();
  int n = Rtime2.size();
  std::vector<int> scoreHCD;
  std::vector<double> rtime1;
  std::vector<double> rtime2;
  std::vector<int> scoreoxonium;
  std::vector<int> mgfi;
  std::vector<int> sizei;
  std::vector<int> sizej;
  std::vector<int> scanindex2;
  std::vector<double> deltamass;
  std::vector<double> neutralloss;
  std::vector< std::string > peptideseq;
  
  NumericVector spectraoxonmium = NumericVector::create(
    163.0606,168.066,186.0166,204.0872,274.0926,292.1032,325.1134,366.14);
  
  NumericVector gapions = NumericVector::create(
   66.0212,73.0259,81.0264,88.0161, 101.5397,132.0423,145.5477,
   146.0579,153.5451,162.0528,176.0321,203.0794,291.0954,307.0903);
  
  for(int i=1;i<m-3;i++)
  {
    if((pepmz1[i]==pepmz1[i+1])&&(pepmz1[i]!=pepmz1[i-1])){
    NumericVector spectraHCD=getspectra(mgfindexexp,mzexp,mgfindex1[i]);
    int scoretempoxo=spectracomp(spectraHCD,spectraoxonmium,0.02);
    if(scoretempoxo>(1)){
    int pos300hcd=mass300pos(spectraHCD);//**mass below 300 are filtered
   
    spectraHCD.erase(spectraHCD.begin(),spectraHCD.begin()+pos300hcd);//**mass below 300 are filtered
    

    double pepmassi=(pepmz1[i]-1.0078)*charge1[i];
    for(int j=0;j<n;j++)
    {
      if((Rtime1[i]-Rtime2[j])<=(-480))
      {
        continue;
      }
      else if((-480<(Rtime1[i]-Rtime2[j]))&&((Rtime1[i]-Rtime2[j])<120))
      {
        std::string seqj=sequence[j];
        double pepmassj=pepmass2[j]-1.0078;
        double tempdeltamass=pepmassi-pepmassj;
        if((tempdeltamass>=202)&(tempdeltamass<=2000))
        {
        NumericVector theospectraHCD=singlespectragen(seqj);
        theospectraHCD=vectorsort(theospectraHCD);  
        int scoretemp=spectracomp(spectraHCD,theospectraHCD,0.02);
        
        if((scoretemp>4)){
        scoreHCD.push_back(scoretemp);
        scoreoxonium.push_back(scoretempoxo);
        deltamass.push_back(tempdeltamass);
        mgfi.push_back(mgfindex1[i]);
        sizei.push_back(spectraHCD.size());
        sizej.push_back(theospectraHCD.size());
        rtime1.push_back(Rtime1[i]);
        rtime2.push_back(Rtime2[j]);
        peptideseq.push_back(seqj);
        scanindex2.push_back(scanindex[j]);
        }
       }
      }
      else continue;
      }
    }
  }
  }
  return Rcpp::DataFrame::create(Named("mgfindex1")=mgfi,Named("sizeHCD")=sizei,Named("sizenative")=sizej, Named("scoreHCD")=scoreHCD,Named("peptideseq")=peptideseq,Named("deltamass")=deltamass,
                                 Named("scoreoxonium")=scoreoxonium, Named("rtimeexp")=rtime1,Named("rtimetheo")=rtime2,Named("scanindex2")=scanindex2);
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
  spectralist1<-spectralist1[intensity>1000,]
  spectralist1<-spectralist1[,head(.SD,200),by=mgfindex]
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
  
  

  spectralist1mgf<-spectralist1$mgfindex
  spectralist1mz<-spectralist1$mz
  spectralist1intensity<-spectralist1$intensity
  isotope<-deisotope(spectralist1mgf,spectralist1mz,spectralist1intensity)
  fragz<-calz(spectralist1mgf,spectralist1mz,spectralist1intensity)
  spectralisttemp<-data.frame(spectralist1,isotope,fragz)
  index<-which(spectralisttemp$isotope==0)
  spectralist1<-spectralisttemp[index,]
  spectralisttemp$mz<-spectralisttemp$mz*spectralisttemp$fragz
  spectralist1<-spectralisttemp[mz<1800,]
  spectralist1<-spectralist1[,1:4]
  mgfindex11<-as.numeric(spectralist1$mgfindex)
  mz1<-as.numeric(spectralist1$mz)
  intensity1<-as.numeric(spectralist1$intensity)

  psms<-read.csv("guoyu-C+F-ID-3-psm.csv")
  psms<-as.data.table(psms)
  indexpsms<-psms[ ,.I[which.max(IDscore)], by=pepsequen]
  indexpsms<-indexpsms$V1
  psms<-psms[indexpsms,.(pepsequen,scanindex,Rtime)]
  psms<-as.data.table(psms)
  setkey(psms,scanindex)
  scanindex2<-psms$scanindex
  sequence<-as.character(psms$pepsequen)
  rtime2<-as.numeric(psms$Rtime)*60
  pepmass2<-calmass(sequence)
  
  
score<-Tspectracomp(mgfindex1,pepmz1,Rtime1,charge1,mgfindex11,mz1,pepmass2,scanindex2,rtime2,sequence)


leng<-length(score$mgfindex1)
scanindex1<-rep(NA,times=leng)
scanindex2<-rep(NA,times=leng)
peptideseq<-rep(NA,times=leng)

for(i in 1:leng)
{
  mgfindex1temp<-score$mgfindex1[i]
  indextemp1<-which(parention1$mgfindex==mgfindex1temp)
  scanindex1[i]<-parention1$scanindex[indextemp1]
  }

score<-data.frame(score,scanindex1)

probscoreHCD<-1-ppois(score$scoreHCD,(score$sizeHCD*score$sizenative/2000))
probscoreHCD<-(-7.2)*log10(probscoreHCD)


newspectra<-data.frame(score,probscoreHCD)
  


leng<-length(newspectra$mgfindex1)
oxoionsDT<-rep(NA,times=leng)
simlarity<-rep(0,times=leng)
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

newspectra<-data.frame(newspectra,oxoionsDT)


write.csv(newspectra,file="score.csv")
  
*/
