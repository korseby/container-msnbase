#!/usr/bin/env Rscript

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)

# Taking the command line arguments
args <- commandArgs(trailingOnly = TRUE)

if(length(args)==0)stop("No file has been specified! Please select a file for performing RT correction!\n")

inputMS2<-NA
inputCamera<-NA
output<-NA
maxSpectra<-NA
minPeaks<-0
minPrecursorMass<-NA
maxPrecursorMass<-NA
precursortppm<-10
fragmentppm<-10
fragmentabs<-0.001
database<-"KEGG"
for(arg in args)
{
  argCase<-strsplit(x = arg,split = "=")[[1]][1]
  value<-strsplit(x = arg,split = "=")[[1]][2]
  
  if(argCase=="inputMS2")
  {
    inputMS2=as.character(value)
  
  }
  if(argCase=="inputCAMERA")
  {
    
    inputCamera=as.character(value)
    
  }
  if(argCase=="maxPrecursorMass")
  {

    maxPrecursorMass=as.numeric(value)

  }
  if(argCase=="minPrecursorMass")
  {

    minPrecursorMass=as.numeric(value)

  }
  if(argCase=="precursorppm")
  {
    
    precursortppm=as.numeric(value)
    
  }
  if(argCase=="fragmentppm")
  {
    
    fragmentppm=as.numeric(value)
    
  }
  
  if(argCase=="fragmentabs")
  {
    
    fragmentabs=as.numeric(value)
    
  }
  if(argCase=="database")
  {
    
    database=as.character(value)
    
  }
  if(argCase=="minPeaks")
  {

    minPeaks=as.numeric(value)

  }
  if(argCase=="maxSpectra")
  {

    maxSpectra=as.numeric(value)

  }
  if(argCase=="output")
  {
    output=as.character(value)
  }
}

if(is.na(inputMS2) | is.na(inputCamera) | is.na(output)) stop("Both input (CAMERA and MS2) and output need to be specified!\n")

settingsObject<-list()
settingsObject[["DatabaseSearchRelativeMassDeviation"]]<-precursortppm
settingsObject[["FragmentPeakMatchAbsoluteMassDeviation"]]<-fragmentabs
settingsObject[["FragmentPeakMatchRelativeMassDeviation"]]<-fragmentppm
settingsObject[["MetFragDatabaseType"]]<-database
settingsObject[["MetFragScoreTypes"]]<-"FragmenterScore"
load(inputMS2)
MappedMS2s<-get(varNameForNextStep)
load(inputCamera)
cameraObject<-get(varNameForNextStep)


source("/usr/local/bin/adductCalculator.r")
source("/usr/local/bin/toMetfragCommand.r")
library(stringr)
toMetfragCommand(mappedMS2 = MappedMS2s$mapped,unmappedMS2 = MappedMS2s$unmapped,
                 cameraObject = cameraObject,searchMultipleChargeAdducts = T,includeUnmapped = F,
                 includeMapped = T,settingsObject = settingsObject,preprocess = F,savePath=output, minPeaks=minPeaks, 
		 maxSpectra=maxSpectra, maxPrecursorMass = maxPrecursorMass, minPrecursorMass = minPrecursorMass)
