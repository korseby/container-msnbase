#!/usr/bin/env Rscript

#adductCalculator(mz=853.33089,charge=2)

adductCalculator<-function(mz=NA,charge=NA,mode="pos",adduct=NA)
{
  
  
  Ion.name=c("M+3H",
             "M+2H+Na",
             "M+H+2Na",
             "M+3Na",
             "M+2H",
             "M+H+NH4",
             "M+H+Na",
             "M+H+K",
             "M+ACN+2H",
             "M+2Na",
             "M+2ACN+2H",
             "M+3ACN+2H",
             "M+H",
             "M+NH4",
             "M+Na",
             "M+CH3OH+H",
             "M+K",
             "M+ACN+H",
             "M+2Na-H",
             "M+IsoProp+H",
             "M+ACN+Na",
             "M+2K-H",
             "M+DMSO+H",
             "M+2ACN+H",
             "M+IsoProp+Na+H",
             "2M+H",
             "2M+NH4",
             "2M+Na",
             "2M+K",
             "2M+ACN+H",
             "2M+ACN+Na")
  
  Ion.mass=c(
    "M/3 + 1.007276",
    "M/3 + 8.334590",
    "M/3 + 15.7661904",
    "M/3 + 22.989218",
    "M/2 + 1.007276",
    "M/2 + 9.520550",
    "M/2 + 11.998247",
    "M/2 + 19.985217",
    "M/2 + 21.520550",
    "M/2 + 22.989218",
    "M/2 + 42.033823",
    "M/2 + 62.547097",
    "M + 1.007276",
    "M + 18.033823",
    "M + 22.989218",
    "M + 33.033489",
    "M + 38.963158",
    "M + 42.033823",
    "M + 44.971160",
    "M + 61.06534",
    "M + 64.015765",
    "M + 76.919040",
    "M + 79.02122",
    "M + 83.060370",
    "M + 84.05511",
    "2M + 1.007276",
    "2M + 18.033823",
    "2M + 22.989218",
    "2M + 38.963158",
    "2M + 42.033823",
    "2M + 64.015765")
  
  Charge=c(
    "3+",
    "3+",
    "3+",
    "3+",
    "2+",
    "2+",
    "2+",
    "2+",
    "2+",
    "2+",
    "2+",
    "2+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+")
  
  adductsFile<-data.frame(Ion.name=Ion.name,Ion.mass=Ion.mass,Charge=Charge)#read.csv("adductList.csv")
  if(is.na(mz))stop("Provide the mz!")
  tmpAdduct<-adductsFile
  if(!is.na(charge) & is.numeric(charge))
  {
    chargeSign<-NA
      if(mode=="pos"){chargeSign<-"+"}else{chargeSign<-"-"}
    tmpAdduct<-tmpAdduct[tmpAdduct[,"Charge"]==paste(charge,chargeSign,sep=""),]
  }
  if(!is.na(adduct))
  {
    
    tmpAdduct<-tmpAdduct[tmpAdduct[,"Ion.name"]%in%adduct,]
  }
  
  if(nrow(tmpAdduct)==0){
    warning("Not found! Reporting all possible adducts")
    tmpAdduct<-adductsFile
    }
  signs<-sapply(str_extract(as.character(tmpAdduct[,"Ion.mass"]),"\\+|\\-"),function(x){x[[1]]})
  signs[signs=="+"]=1
  signs[signs=="-"]=-1
  result<-(mz*as.numeric(tmpAdduct[,"Charge"]))-
    (as.numeric(signs)*
       as.numeric(sapply(strsplit(as.character(tmpAdduct[,"Ion.mass"]),"\\+|\\-"),function(x){x[[2]]})))
  
return(data.frame(correctedMS=result,adductName=tmpAdduct[,"Ion.name"]))
  
}
