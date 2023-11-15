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

std::vector<double>  singlespectragen(std::string pepsequen1){
  
  std::vector< double > mz;
  std::map<char,double> AAmass;
  AAmass.insert ( std::pair<char,double>('A',71.03711) );
  AAmass.insert ( std::pair<char,double>('R',156.10111) );
  AAmass.insert ( std::pair<char,double>('N',114.04293) );
  AAmass.insert ( std::pair<char,double>('D',115.02694) );
  AAmass.insert ( std::pair<char,double>('C',160.0338) );
  AAmass.insert ( std::pair<char,double>('E',129.04259) );
  AAmass.insert ( std::pair<char,double>('Q',128.05858) );
  AAmass.insert ( std::pair<char,double>('G',57.02146) );
  AAmass.insert ( std::pair<char,double>('H',137.05891) );
  AAmass.insert ( std::pair<char,double>('I',113.08406) );
  AAmass.insert ( std::pair<char,double>('L',113.08406) );
  AAmass.insert ( std::pair<char,double>('K',128.09496) );
  AAmass.insert ( std::pair<char,double>('M',131.04049) );
  AAmass.insert ( std::pair<char,double>('F',147.06841) );
  AAmass.insert ( std::pair<char,double>('P',97.05276) );
  AAmass.insert ( std::pair<char,double>('S',87.03203) );
  AAmass.insert ( std::pair<char,double>('T',101.04768) );
  AAmass.insert ( std::pair<char,double>('W',186.07931) );
  AAmass.insert ( std::pair<char,double>('Y',163.06333) );
  AAmass.insert ( std::pair<char,double>('V',99.06841) );
  
  
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
  return mz;
}






// [[Rcpp::export]]

NumericVector  calmass(std::vector< std::string > seq1){
  
  
  int m=seq1.size();
  
  std::vector< double > pepmass;
 
  
  std::map<char,double> AAmass;
  AAmass.insert ( std::pair<char,double>('A',71.03711) );
  AAmass.insert ( std::pair<char,double>('R',156.10111) );
  AAmass.insert ( std::pair<char,double>('N',114.04293) );
  AAmass.insert ( std::pair<char,double>('D',115.02694) );
  AAmass.insert ( std::pair<char,double>('C',160.0338) );
  AAmass.insert ( std::pair<char,double>('E',129.04259) );
  AAmass.insert ( std::pair<char,double>('Q',128.05858) );
  AAmass.insert ( std::pair<char,double>('G',57.02146) );
  AAmass.insert ( std::pair<char,double>('H',137.05891) );
  AAmass.insert ( std::pair<char,double>('I',113.08406) );
  AAmass.insert ( std::pair<char,double>('L',113.08406) );
  AAmass.insert ( std::pair<char,double>('K',128.09496) );
  AAmass.insert ( std::pair<char,double>('M',131.04049) );
  AAmass.insert ( std::pair<char,double>('F',147.06841) );
  AAmass.insert ( std::pair<char,double>('P',97.05276) );
  AAmass.insert ( std::pair<char,double>('S',87.03203) );
  AAmass.insert ( std::pair<char,double>('T',101.04768) );
  AAmass.insert ( std::pair<char,double>('W',186.07931) );
  AAmass.insert ( std::pair<char,double>('Y',163.06333) );
  AAmass.insert ( std::pair<char,double>('V',99.06841) );
  
  
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

int mass300pos(std::vector<double> x){
  int m=x.size();
  int i;
  for(i=0;i<m;i++){if(x[i]>200) break;}
  
  
  
  return(i);
  
}

// [[Rcpp::export]]

std::vector<double> ADD(std::vector<double> x, double delta){
  int m=x.size();

  std::vector<double> mz;
  for(int i=0;i<m;i++)
    {mz.push_back(x[i]+delta);}
  
  return mz;
  
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
      if((y[i]>350)&&((fabs(y[j]-y[i]-1.0034) <0.01) |(fabs(y[j]-y[i]-0.50017)<0.01)|(fabs(y[j]-y[i]-0.3344)<0.01)|(fabs(y[j]-y[i]-0.2509)<0.006)|(fabs(y[j]-y[i]-0.20068)<0.006))&&(x[i]==x[j])&&(z[i]>2000)&&(z[j]/z[i]>y[i]/8000)&&(z[j]/z[i]<y[i]/300))
        iso[j]=1;  
    }
  }
  
  return wrap(iso); 
  
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
      if((y[i]>250)&&((fabs(y[j]-y[i]-0.50017)<0.01))&&(x[i]==x[j])&&(z[i]>2000)&&(z[j]/z[i]>y[i]/8000)&&(z[j]/z[i]<y[i]/300))
        iso[i]=2;  
    }
  }
  
  
  for(int i=0;i<m;i++)
  {
    int minsize=std::min(i+5,m-1);
    for(int j=i+1;j<minsize;j++)
    {
      if((y[i]>250)&&((fabs(y[j]-y[i]-0.3344) <0.01) )&&(x[i]==x[j])&&(z[i]>2000)&&(z[j]/z[i]>y[i]/8000)&&(z[j]/z[i]<y[i]/300))
        iso[i]=3;  
    }
  }
  
  for(int i=0;i<m;i++)
  {
    int minsize=std::min(i+5,m-1);
    for(int j=i+1;j<minsize;j++)
    {
      if((y[i]>250)&&((fabs(y[j]-y[i]-0.25085) <0.01) )&&(x[i]==x[j])&&(z[i]>2000)&&(z[j]/z[i]>y[i]/8000)&&(z[j]/z[i]<y[i]/300))
        iso[i]=4;  
    }
  }
  for(int i=0;i<m;i++)
  {
    int minsize=std::min(i+5,m-1);
    for(int j=i+1;j<minsize;j++)
    {
      if((y[i]>250)&&((fabs(y[j]-y[i]-0.20068) <0.01) )&&(x[i]==x[j])&&(z[i]>2000)&&(z[j]/z[i]>y[i]/8000)&&(z[j]/z[i]<y[i]/300))
        iso[i]=5;  
    }
  }
  
  return wrap(iso);
  
}


// [[Rcpp::export]]

NumericVector posgetspectra(NumericVector x, NumericVector y,int t,int pos){
  
  NumericVector z;
  int m=x.size();
  int tempos=std::max(0,pos-2);
  
  for(int i=tempos;i<m;i++)
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

std::vector<double> posgetintenspectra(NumericVector x, NumericVector y,NumericVector z,int t,int pos){
  
  std::vector<double> mz1;
  std::vector<double> mz;
  std::vector<double> inten;
  int m=x.size();
  int tempos=std::max(0,pos-2);
  
  for(int i=tempos;i<m;i++)
  {
    if ( x[i]<t) 
    {
      continue;
    }
    
    else if ( x[i]==t) 
    {
      mz1.push_back(y[i]);
      inten.push_back(z[i]);
    }
    else
      break;
  }
  
  for(int i=0;i<mz1.size();i++){
    if(inten[i]>2000){
      mz.push_back(mz1[i]);
    }
  }
  return mz;
  
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

std::vector<double> getspectra(NumericVector x, NumericVector y,int t){
  
  std::vector<double> z;
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
     
  return z;
  
}

// [[Rcpp::export]]

std::vector<double> getintenspectra(NumericVector x, NumericVector y,NumericVector z,int t){
  
  std::vector<double> mz1;
  std::vector<double> mz;
  std::vector<double> inten;
  int m=x.size();
  
  
  for(int i=0;i<m;i++)
  {
    if ( x[i]<t) 
    {
      continue;
    }
    
    else if ( x[i]==t) 
    {
      mz1.push_back(y[i]);
      inten.push_back(z[i]);
    }
    else
      break;
  }
  
  for(int i=0;i<mz1.size();i++){
    if(inten[i]>2000){
      mz.push_back(mz1[i]);
    }
  }
  return mz;
  
  
}

// [[Rcpp::export]]

int spectracomp(std::vector<double> x, NumericVector y,double delta){
  
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

int spectracomp2(std::vector<double> x, NumericVector y,double delta,double eppm){
  
  int m = x.size();
  int n = y.size(); 
  
  int i=0;
  int j=0;
  std::vector<double> deltamass;
  
  
  while(i<m && j<n)
  {
    
    if ((x[i]-y[j])<delta&&(x[i]-y[j])>(-delta))
    {
      double xtemp=x[i];
      double ytemp=y[j];
      double error=fabs(xtemp-ytemp);
      error=error/xtemp;
      deltamass.push_back(error);
      i=i+1;
      j=j+1;
    }
    else if(x[i]<=y[j]-delta) i=i+1;
    else if(x[i]>=y[j]+delta) j=j+1;
  }
  
  int t=0;
  int indexsize=deltamass.size();
  for(int w=0;w<indexsize;w++){
    
    if((deltamass[w]*1e6)<eppm){
      t++;}
    else continue;
  }
  
  return t;
}




// [[Rcpp::export]]

int spectracomp3(std::vector<double> x, std::vector<double> y,double delta){
  
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

int spectracomp4(std::vector<double> x, std::vector<double> y,double delta,double eppm){
  
  int m = x.size();
  int n = y.size(); 
  
  int i=0;
  int j=0;
  std::vector<double> deltamass;
  
  
  while(i<m && j<n)
  {
    
    if ((x[i]-y[j])<delta&&(x[i]-y[j])>(-delta))
    {
      double xtemp=x[i];
      double ytemp=y[j];
      double error=fabs(xtemp-ytemp);
      error=error/xtemp;
      deltamass.push_back(error);
      i=i+1;
      j=j+1;
    }
    else if(x[i]<=y[j]-delta) i=i+1;
    else if(x[i]>=y[j]+delta) j=j+1;
  }
  
  int t=0;
  int indexsize=deltamass.size();
  for(int w=0;w<indexsize;w++){
    
    if((deltamass[w]*1e6)<eppm){
      t++;}
    else continue;
  }
  
  return t;
}




// [[Rcpp::export]]

DataFrame Tspectracomp(NumericVector mgfindex1,NumericVector pepmz1,IntegerVector Rtime1,NumericVector charge1,
                        NumericVector mgfindexexp,NumericVector mzexp,NumericVector pepmass2,IntegerVector scanindex,
                        NumericVector Rtime2,std::vector< std::string > sequence) {
  
  
 
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
  std::vector<double> pepmz;
  std::vector<int> charge;
  std::vector< std::string > peptideseq;

  
  NumericVector spectraoxonmium = NumericVector::create(
    126.0550,138.0550,145.0500,163.0606,168.066,186.0166,
    204.0872,274.0926,292.1032,325.1134,366.144);
  
  NumericVector gapions = NumericVector::create(
   66.0212,73.0259,81.0264,88.0161, 101.5397,132.0423,145.5477,
   146.0579,153.5451,162.0528,176.0321,203.0794,291.0954,307.0903);
  
  for(int i=1;i<m-3;i++)
  {
    std::vector<double> expspectraHCD1;
    std::vector<double> expspectraHCD2;
    
    
    if((pepmz1[i]==pepmz1[i+1])&&(pepmz1[i]!=pepmz1[i-1])){
      expspectraHCD1=getspectra(mgfindexexp,mzexp,mgfindex1[i]);
      expspectraHCD2=getspectra(mgfindexexp,mzexp,mgfindex1[i]+1);
      }
    
    else if((pepmz1[i]!=pepmz1[i+1])&&(pepmz1[i]==pepmz1[i-1])){
      expspectraHCD1=getspectra(mgfindexexp,mzexp,mgfindex1[i]+1);
      expspectraHCD2=getspectra(mgfindexexp,mzexp,mgfindex1[i]);
      
    }
    
    
    
    
    std::vector<double> spectraHCD=getspectra(mgfindexexp,mzexp,mgfindex1[i]);
    int scoretempoxo=spectracomp2(expspectraHCD1,spectraoxonmium,0.02,20);
    if(scoretempoxo>(1)){
    int pos300hcd=mass300pos(expspectraHCD2);//**mass below 300 are filtered
   
    spectraHCD.erase(expspectraHCD2.begin(),expspectraHCD2.begin()+pos300hcd);//**mass below 300 are filtered
    

    double pepmassi=(pepmz1[i]-1.0078)*charge1[i];
    for(int j=0;j<n;j++)
    {
      if((Rtime1[i]-Rtime2[j])<=(-480))
      {
        continue;
      }
      else if((-480<(Rtime1[i]-Rtime2[j]))&&((Rtime1[i]-Rtime2[j])<60))
      {
        std::string seqj=sequence[j];
        double pepmassj=pepmass2[j]-1.0078;
        double tempdeltamass=pepmassi-pepmassj;
        if((tempdeltamass>=(120))&(tempdeltamass<=2000))
        {
          std::vector<double> theospectraHCD=singlespectragen(seqj);
       
        int scoretemp=spectracomp4(expspectraHCD2,theospectraHCD,0.06,20);
        
        if((scoretemp>4)){
        scoreHCD.push_back(scoretemp);
        scoreoxonium.push_back(scoretempoxo);
        deltamass.push_back(tempdeltamass);
        mgfi.push_back(mgfindex1[i]);
        sizei.push_back(expspectraHCD2.size());
        sizej.push_back(theospectraHCD.size());
        rtime1.push_back(Rtime1[i]);
        rtime2.push_back(Rtime2[j]);
        peptideseq.push_back(seqj);
        scanindex2.push_back(scanindex[j]);
        pepmz.push_back(pepmz1[i]);
        charge.push_back(charge1[i]);
       
        }
       }
     
      else continue;
      }
    }
  }
  }
  return Rcpp::DataFrame::create(Named("Gmgfindex1")=mgfi,Named("sizeexp")=sizei,Named("sizeHCD")=sizej, Named("scoreHCD")=scoreHCD,Named("PeptideSequence")=peptideseq,Named("Glycanmass")=deltamass,
                                 Named("scoreoxonium")=scoreoxonium, Named("GlycoRtime")=rtime1,Named("IDritme")=rtime2,Named("GLycoscan")=scanindex2,Named("pepmz")=pepmz,Named("GLycocharge")=charge);
}



// [[Rcpp::export]]

DataFrame similaritysearch(NumericVector gmgfindex1,NumericVector gRtime1,NumericVector gscore1, std::vector<std::string> sequence,NumericVector gmass,
                         NumericVector mgfindex1,NumericVector pepmz1,NumericVector Rtime1,NumericVector charge1,IntegerVector mgfpos,
                         NumericVector mgfindex11,NumericVector mz1,NumericVector intensity1) {
 

  
  

  int m = gmgfindex1.size();
  int n= mgfindex1.size();
  std::vector<int> score;
  std::vector<int> scoreoxonium;
  std::vector<double> gscore;
  std::vector<int> mgfi;
  std::vector<int> mgfj;
  std::vector<int> sizei;
  std::vector<int> sizej;
  std::vector<double> deltamass;
  std::vector< std::string > peptideseq;
  std::vector<double> pepmz;
  std::vector<int> charge;
  
  NumericVector spectraoxonmium = NumericVector::create(
    126.0550,138.0550,145.0500,163.0606,168.066,186.0166,
    204.0872,274.0926,292.1032,325.1134,366.14);
  
  
  
  for(int i=0;i<m;i++)
  {
    std::vector<double> spectra1=getintenspectra(mgfindex11,mz1,intensity1,gmgfindex1[i]);
   
    double pepmassi=(pepmz1[i]-1.0078)*charge1[i];
    int pos3001=mass300pos(spectra1);//**mass below 300 are filtered
    spectra1.erase(spectra1.begin(),spectra1.begin()+pos3001);//**mass below 300 are filtered
    
    for(int j=1;j<n;j++)
    {
     if((-240<(Rtime1[j]-gRtime1[i]))&&((Rtime1[j]-gRtime1[i])<240))
      {
        std::string seqi=sequence[i];
        std::vector<double> spectra2;
        
        if((pepmz1[j]==pepmz1[j+1])&&(pepmz1[j]!=pepmz1[j-1])){
          spectra2=posgetintenspectra(mgfindex11,mz1,intensity1,mgfindex1[j],mgfpos[j]);
        }
        
    
        double deltamassi=gmass[j]-pepmassi;
        int scoretempoxo=spectracomp2(spectra2,spectraoxonmium,0.02,20);
        int pos3002=mass300pos(spectra2);//**mass below 300 are filtered
        spectra2.erase(spectra2.begin(),spectra2.begin()+pos3002);//**mass below 300 are filtered
        std::vector<double> spectra2delta=ADD(spectra2,deltamassi);
        int scoretemp=spectracomp3(spectra1,spectra2,0.05);
          //+spectracomp3(spectra1,spectra2delta,0.06);
        
         if(scoretemp>4&&scoretempoxo>1){
          score.push_back(scoretemp);
          scoreoxonium.push_back(scoretempoxo);
          gscore.push_back(gscore1[i]);
          mgfi.push_back(gmgfindex1[i]);
          mgfj.push_back(mgfindex1[j]);
          sizei.push_back(spectra1.size());
          sizej.push_back(spectra2.size());
          pepmz.push_back(pepmz1[j]);
          charge.push_back(charge1[j]);
          peptideseq.push_back(seqi);
          deltamass.push_back(deltamassi);
        }
      }
      else continue;
    }
    }
  return Rcpp::DataFrame::create(Named("G1mgfindex")=mgfi,Named("size1")=sizei,Named("size2")=sizej,
                                    Named("scoreoxonium")=scoreoxonium,Named("score")=score,Named("pepmz")=pepmz,Named("charge")=charge,
                                       Named("G2mgfindex")=mgfj,Named("PeptideSequence")=peptideseq,Named("Gscore")=gscore,Named("deltamass")=deltamass);
}

/*** R
setwd("D:\\Rfiles\\Glycan")
library(data.table)
  library(Rcpp)
  library(inline)
  library(ggplot2)
  
  
  
mgfread<-function(df)
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
}# readmgf files and convert the mgf files to a list of mz values including charge and scanindex et al.


fragmentprocessor<-function(spectralist){
  spectralist1<-as.data.table(spectralist)
  spectralist1$mz<-round(as.numeric(spectralist1$mz),4)
  spectralist1<-spectralist1[,.SD[1],by=.(mgfindex,mz)]
  setkey(spectralist1,mgfindex,intensity)
  spectralist1<-spectralist1[intensity>500,]
  spectralist1<-spectralist1[,head(.SD,300),by=mgfindex]
  mgfindex<-as.numeric(spectralist1$mgfindex)
  mz<-as.numeric(spectralist1$mz)
  spectralist1<-spectralist1[order(mgfindex,mz),]
  spectralist1mgf<-as.numeric(spectralist1$mgfindex)
  spectralist1mz<-as.numeric(spectralist1$mz)
  spectralist1intensity<-as.numeric(spectralist1$intensity)
  isotope<-deisotope(spectralist1mgf,spectralist1mz,spectralist1intensity)
  fragz<-calz(spectralist1mgf,spectralist1mz,spectralist1intensity)
  spectralisttemp<-data.frame(spectralist1,isotope,fragz)
  spectralisttemp$mz<-(spectralisttemp$mz-1.0078)*spectralisttemp$fragz+1.0078
  index<-which(spectralisttemp$isotope==0)
  spectralist1<-spectralisttemp[index,]
  spectralist1<-spectralist1[,1:6]
  spectralist1<-as.data.table(spectralist1)
  spectralist1<-spectralist1[,.SD[1],by=.(mgfindex,mz)]
  setkey(spectralist1,mgfindex,mz)
  mgfindex11<-as.numeric(spectralist1$mgfindex)
  mz1<-as.numeric(spectralist1$mz)
  intensity1<-as.numeric(spectralist1$intensity)
  spectrapos<-spectralist1[,.N,by=mgfindex]
  mgfleng<-spectrapos$N
  mgfpos<-cumsum(mgfleng)
  mgfpos<-c(0,mgfpos)
  leng<-length(mgfpos)
  mgfpos<-mgfpos[-leng]
  mgfpos<-as.integer(mgfpos)
  mgfindex<-spectrapos$mgfindex
  return (list("mgfindex11"=mgfindex11,"mz1"=mz1,"intensity1"=intensity1,"mgfpos"=mgfpos,"mgfleng"=mgfleng,"mgfindex"=mgfindex))
}  #Generate fragmentation mz and intensity information for each scan






  #sax<-fread("jiaoda-IGA-DYC-4.mgf",sep="=", fill=TRUE,blank.lines.skip=TRUE,data.table=TRUE)
  #sax<-data.frame(sax)
  #spectra<-mgfread(sax)
  #parention1<-spectra$parention
  #spectralist1<-spectra$fragmention
  #write.csv(spectralist1,file="spectralist1.csv")
  #write.csv(parention1,file="parention1.csv")
  #rm(spectra,spectralist1,parention1)

  parention1<-read.csv("parention1.csv")
  spectralist1<-fread("spectralist1.csv")
  spectralist<-fragmentprocessor(spectralist1)
  mgfindex11<- spectralist$mgfindex11
  mz1<-spectralist$mz1
  intensity1<-spectralist$intensity1
  mgfpos<-spectralist$mgfpos
  mgfleng<-spectralist$mgfleng
  mgfindex<-spectralist$mgfindex
  parention2<-data.frame(mgfindex,mgfpos,mgfleng)
  parention1<-data.table(parention1)
  setkey(parention1,mgfindex)
  parention2<-data.table(parention2)
  setkey(parention2,mgfindex)
  parention<-merge(parention1,parention2,by="mgfindex",all=TRUE)
  write.csv(parention,file="parention.csv")
  rm(parention)
  parention<-read.csv("parention.csv")
  parention<-as.data.table(parention)
  setkey(parention,mgfindex,scanindex)
  mgfindex1<-as.integer(parention$mgfindex)
  pepmz1<-parention$pepmass
  Rtime1<-parention$Rtime
  charge1<-as.integer(parention$charge)
  scanindex1<-as.integer(parention$scanindex)
  pepmass1<-(pepmz1-1.0078)*charge1+1.0078
  mgfpos<-as.integer(parention$mgfpos)
  mgfleng<-as.integer(parention$mgfleng)
  index<-which(is.na(mgfpos))
  mgfpos[index]<-length(mz1)-1
  mgfleng[index]<-0



  psms<-read.csv("jiaoda-IGA-DYC-4-PSMS.csv")
  psms<-as.data.table(psms)
  psms<-psms[,mRtime:=median(Rtime),by=AnnotatedSequence]
  indexpsms<-psms[ ,.I[which.max(IDScore)], by=AnnotatedSequence]
  indexpsms<-indexpsms$V1
  psms<-psms[indexpsms,.(ProteinAccessions,AnnotatedSequence,FirstScan,mRtime)]
  psms<-as.data.table(psms)
  setkey(psms,FirstScan)
  scanindex2<-psms$FirstScan
  sequence<-as.character(psms$AnnotatedSequence)
  rtime2<-as.numeric(psms$mRtime)*60
  pepmass2<-calmass(sequence)
  
  
  t1<-Sys.time()
  
 score<-Tspectracomp(mgfindex1,pepmz1,Rtime1,charge1,mgfindex11,mz1,pepmass2,scanindex2,rtime2,sequence)


 leng<-length(score$Gmgfindex1)
 IDscan<-rep(NA,times=leng)

  peptideseq<-rep(NA,times=leng)

  for(i in 1:leng)
{
 mgfindex1temp<-score$Gmgfindex1[i]
 IDscan[i]<-parention1$scanindex[mgfindex1temp]
 }

score<-data.frame(score,IDscan)

probscoreHCD<-1-ppois(score$scoreHCD,(score$sizeHCD*score$sizeexp/2000))
GlycoScore<-round((-7.2)*log10(probscoreHCD),2)


newspectra<-data.frame(score,GlycoScore)
index<-which(newspectra$GlycoScore>10)
newspectra<-newspectra[index,]

newspectra$PeptideSequence<-as.character(newspectra$PeptideSequence)
glycanmass<-newspectra$GlycanMass+calmass(newspectra$PeptideSequence)-1.0078

newspectra2<-similaritysearch(newspectra$Gmgfindex1,newspectra$GlycoRtime,newspectra$GlycoScore,
                              newspectra$PeptideSequence,glycanmass,
                                  mgfindex1,pepmz1,Rtime1,charge1,mgfpos,mgfindex11,mz1,as.integer(intensity1))

t2<-Sys.time()
print(t2-t1)
leng<-length(newspectra2$G1mgfindex)
Protein_Name<-rep(NA,times=leng)
G1scanindex<-rep(NA,times=leng)
G2scanindex<-rep(NA,times=leng)
IDRtime<-rep(NA,times=leng)
G2Rtime<-rep(NA,times=leng)

for(i in 1:leng)
{
  seqtemp<-as.character(newspectra2$PeptideSequence[i])
  index<-which(psms$AnnotatedSequence==seqtemp)
  index1<-index[1]
  Protein_Name[i]<-as.character(psms$ProteinAccessions[index1])
  
  
  mgfindex1i<-newspectra2$G1mgfindex[i]
  mgfindex2i<-newspectra2$G2mgfindex[i]
  G1scanindex[i]<-parention$scanindex[mgfindex1i]
  G2scanindex[i]<-parention$scanindex[mgfindex2i]
  IDRtime[i]<-median(psms$mRtime[index])*60
  G2Rtime[i]<-parention$Rtime[mgfindex2i]
 }

newspectra2<-data.frame(Protein_Name,G1scanindex,G2scanindex,IDRtime,G2Rtime,newspectra2)
simscore<-1-ppois(newspectra2$score,(newspectra2$size1*newspectra2$size2/2000))
ClusterScore<-round((-7.2)*log10(simscore),2)
Glymass<-(newspectra2$pepmz-1.0078)*newspectra2$charge-calmass(as.character(newspectra2$PeptideSequence))+1.0078
newscore<-data.frame(newspectra2,ClusterScore,Glymass)
index<-which(newscore$Glymass>120)
newscore<-newscore[index,]



write.csv(newscore,file="jiaoda-IGA-DYC-4.csv")
write.csv(newspectra,file="jiaoda-IGA-DYC-4-C1st.csv")
*/
