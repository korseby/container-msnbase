#!/usr/bin/env Rscript
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)

# Taking the command line arguments
args <- commandArgs(trailingOnly = TRUE)

if(length(args)==0)stop("No file has been specified! Please select a file for performing peak picking!\n")
require(MSnbase)
RawFiles<-NA
precursorShift<-0
output<-NA
for(arg in args)
{
  argCase<-strsplit(x = arg,split = "=")[[1]][1]
  value<-strsplit(x = arg,split = "=")[[1]][2]
  if(argCase=="input")
  {
    RawFiles=as.character(value)
  }
  if(argCase=="precursorShift")
  {
    precursorShift=as.numeric(value)
  }
  if(argCase=="output")
  {
    output=as.character(value)
  }
  
}
if(is.na(RawFiles) | is.na(output)) stop("Both input and output need to be specified!\n")
MS2RawFile<-readMSData(RawFiles, msLevel = 2, verbose = FALSE)
sapply(1:length(MS2RawFile), function(x) MS2RawFile[[x]]@precursorMz<<-MS2RawFile[[x]]@precursorMz + precursorShift)

preprocessingStepsMS2<-c("MS2RawFile")
varNameForNextStep<-as.character("MS2RawFile")
save(list = c("MS2RawFile","preprocessingStepsMS2","varNameForNextStep"),file = output)
