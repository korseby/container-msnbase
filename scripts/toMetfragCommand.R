#!/usr/bin/env Rscript

parameterToCommand<-function(param,outputName="")
{
  param$MetFragPeakListReader<-"de.ipbhalle.metfraglib.peaklistreader.FilteredStringTandemMassPeakListReader"
  param$MetFragCandidateWriter<-"CSV"
  param$PeakListString<-paste(apply(param$PeakList,1,paste,collapse="_"),collapse = ";")
  param$SampleName<-outputName
  param[[which(names(param)=="PeakList")]]<-NULL

 param$MetFragDatabaseType
 toOutput<-""
 for(i in 1:length(param))
 {
   if(i==1) {toOutput<- paste(names(param)[i],"=",param[[i]],sep="")
   } else {toOutput<- paste(toOutput," ",names(param)[i],"=",param[[i]],sep="")}
 }
 ### the output is rt_mz_randomNumber. the random number is to prevent overwritting of ms2 withs the same precursor mz
 cat(toOutput,file = outputName)
}
require(CAMERA)
require(stringr)
toMetfragCommand<-function(mappedMS2=NA,
                           unmappedMS2=NA,
                           cameraObject=NA,
                           searchMultipleChargeAdducts=F,
                           includeUnmapped=T,includeMapped=T,
                           settingsObject=list(),preprocess=NA,savePath="",minPeaks=0,maxSpectra=NA,
			   maxPrecursorMass = NA, minPrecursorMass = NA)
{
  print(paste("lens",length(mappedMS2),length(unmappedMS2))
  peakList<-getPeaklist(cameraObject)
  numberSpectraWritten <- 0
  if(includeMapped==T)
  {
    searchChargeFlag<-F
    searchAdductsFlag<-F
    massTraceNames<-names(mappedMS2)
    metFragResult<-c()
    for(x in massTraceNames)
    {
      seachAdducts<-NA
      seachCharge<-NA       
      searchChargeFlag<-F       
      if(peakList[as.numeric(x),"adduct"]=="" & peakList[as.numeric(x),"isotopes"]=="")
      {
        adduct<-NA
        neutralMASS<-peakList[as.numeric(x),"mz"]
        searchChargeFlag<-T
        searchAdductsFlag<-T
        seachAdducts<-adduct
        
      }else if(peakList[as.numeric(x),"adduct"]!="")
      {
        if(str_count(peakList[as.numeric(x),"adduct"], "]")==1)
        {
          adduct<-gsub(" ","",str_extract(peakList[as.numeric(x),"adduct"], "\\[.*\\]*. "))
          neutralMASS<-as.numeric(str_extract(peakList[as.numeric(x),"adduct"], " .*"))
          searchChargeFlag<-F
        }else 
        {
          adduct<-""
          neutralMASS<-peakList[as.numeric(x),"mz"]
          adduct<-unique(sapply(strsplit(str_extract_all(peakList[as.numeric(x),"adduct"], "\\[.*?\\]")[[1]]," "),
                                function(x){x[[1]]}))
          searchChargeFlag<-T
          seachAdducts<-adduct
        }
      }else if(peakList[as.numeric(x),"isotopes"]!="")
      {
        
        isotopID<- str_extract(peakList[as.numeric(x),"isotopes"], "\\[\\d*\\]")
        
        monoIsotopic<-grepl("[M]",peakList[ grepl(isotopID,peakList[,"isotopes"],fixed=T),"isotopes"],fixed=T)
        adduct<- gsub(" ","",
                      str_extract(peakList[ grepl(isotopID,peakList[,"isotopes"],fixed=T),][monoIsotopic,"adduct"],"\\[.*\\]*. "))
        
        tmpMASS<-as.numeric(str_extract(peakList[ grepl(isotopID,peakList[,"isotopes"],fixed=T),][monoIsotopic,"adduct"], " .*"))
        if(adduct!="" & !is.na(adduct))
        {
          tmpMASS<-as.numeric(str_extract(peakList[ grepl(isotopID,peakList[,"isotopes"],fixed=T),][monoIsotopic,"adduct"], " .*"))
          neutralMASS<-tmpMASS			 
        }else
        { 
          isoTMP<-gsub("\\[.*\\]","",peakList[as.numeric(x),"isotopes"])
          charges<-str_extract(isoTMP,"\\d+")
          
          if(is.na(charges))
          {
            adduct<-NA
            neutralMASS<-peakList[as.numeric(x),"mz"]
            searchChargeFlag<-T
            seachCharge<-1
          }else
          {
            neutralMASS<-peakList[as.numeric(x),"mz"]
            searchChargeFlag<-T
            seachCharge<-charges
          }
          
        }
        
      }
      
      mappedMS2TMP<-NA
      if(class(mappedMS2[[x]])=="list")
      {
        mappedMS2TMP<-mappedMS2[[x]]
      }
      else if(class(mappedMS2[[x]])=="Spectrum2")
      {
        mappedMS2TMP<-list(mappedMS2[[x]])
      }
      for(MSMS in mappedMS2TMP)
      {
        if(preprocess==T) 
        {
          
          MSMS@centroided<-F
          MSMS@polarity<-as.integer(1)
          MSMS@smoothed<-F
          MSMS<-MSnbase::pickPeaks(MSMS)
          
        }
        MS2<-as.matrix(cbind(MSMS@mz,MSMS@intensity))
        # if number MS/MS peaks is too low
	print(MSMS@mz)
        print(length(MSMS@mz)
	if(length(MSMS@mz) == 0 || dim(MS2)[1] < minPeaks) { next }
        if(searchChargeFlag==F)
        {
          settingsObject[["NeutralPrecursorMass"]]<-neutralMASS
          settingsObject[["PeakList"]]<-MS2
          fileName<-""
          fileName<-paste(as.character(MSMS@rt),"_",as.character(round(neutralMASS,4)),"_",as.character(runif(1)),".txt",sep="")
         if(savePath!="")
          fileName<-paste(savePath,"/",as.character(MSMS@rt),"_",as.character(round(neutralMASS,4)),"_",as.character(runif(1)),".txt",sep="")
         if(!is.na(maxPrecursorMass) & maxPrecursorMass < neutralMASS) { next }
         if(!is.na(minPrecursorMass) & minPrecursorMass > neutralMASS) { next }
	 if(is.na(maxSpectra) || maxSpectra > numberSpectraWritten) {
         	parameterToCommand(settingsObject,fileName)
	 	numberSpectraWritten<-numberSpectraWritten+1
	 }
        }else if(searchChargeFlag==T & searchMultipleChargeAdducts==T)
        {
          
          allChargesHits<-list()
          allAdductForSearch<-adductCalculator(mz = neutralMASS,charge = seachCharge,
                                               adduct = gsub("\\[|\\]","",seachAdducts),mode = "pos")
          for(k in nrow(allAdductForSearch))
          {
            mass <- allAdductForSearch[k,"correctedMS"]
            settingsObject[["NeutralPrecursorMass"]]<-mass
            settingsObject[["PeakList"]]<-MS2
            fileName<-""
            fileName<-paste(as.character(MSMS@rt),"_",as.character(round(mass,4)),"_",as.character(runif(1)),".txt",sep="")
            if(savePath!="")
              fileName<-paste(savePath,"/",as.character(MSMS@rt),"_",as.character(round(mass,4)),"_",as.character(runif(1)),".txt",sep="")
            if(!is.na(maxPrecursorMass) & maxPrecursorMass < mass) { next }
            if(!is.na(minPrecursorMass) & minPrecursorMass > mass) { next }
            if(is.na(maxSpectra) || maxSpectra > numberSpectraWritten) {
	    	parameterToCommand(settingsObject,fileName)
	    	numberSpectraWritten<-numberSpectraWritten+1
	    }
          }
        }
        
      }
      
    }
  }
  
  if(includeUnmapped)
  {
    
    for(p in 1:length(unmappedMS2))
    {
      MSMS<-unmappedMS2[[p]]
      if(preprocess==T) 
      {
        
        MSMS@centroided<-F
        MSMS@polarity<-as.integer(1)
        MSMS@smoothed<-F
        MSMS<-MSnbase::pickPeaks(MSMS)
        
      }
      neutralMASS<-MSMS@precursorMz
      MS2<-as.matrix(cbind(MSMS@mz,MSMS@intensity))
      print(MSMS@mz)
      print(length(MSMS@mz)
      if(length(MSMS@mz) == 0 || dim(MS2)[1] < minPeaks) { next }
      if(searchMultipleChargeAdducts==F)
      {
        settingsObject[["NeutralPrecursorMass"]]<-neutralMASS
        settingsObject[["PeakList"]]<-MS2
        fileName<-""
        fileName<-paste(as.character(MSMS@rt),"_",as.character(round(neutralMASS,4)),"_",as.character(runif(1)),".txt",sep="")
        if(savePath!="")
          fileName<-paste(savePath,"/",as.character(MSMS@rt),"_",as.character(round(neutralMASS,4)),"_",as.character(runif(1)),".txt",sep="")
        if(!is.na(maxPrecursorMass) & maxPrecursorMass < neutralMASS) { next }
        if(!is.na(minPrecursorMass) & minPrecursorMass > neutralMASS) { next }
	if(is.na(maxSpectra) || maxSpectra > numberSpectraWritten) {
		parameterToCommand(settingsObject,fileName)
        	numberSpectraWritten<-numberSpectraWritten+1 
	}
      }else if(searchMultipleChargeAdducts==T)
      {
        allChargesHits<-list()
        allAdductForSearch<-adductCalculator(mz = neutralMASS,charge = NA,
                                             adduct = NA,mode = "pos")
        for(k in 1:nrow(allAdductForSearch))
        {
          mass <- allAdductForSearch[k,"correctedMS"]
          settingsObject[["NeutralPrecursorMass"]]<-mass
          settingsObject[["PeakList"]]<-MS2
          fileName<-""
          fileName<-paste(as.character(MSMS@rt),"_",as.character(round(mass,4)),"_",as.character(runif(1)),".txt",sep="")
          if(savePath!="")
            fileName<-paste(savePath,"/",as.character(MSMS@rt),"_",as.character(round(mass,4)),"_",as.character(runif(1)),".txt",sep="")
          if(!is.na(maxPrecursorMass) & maxPrecursorMass < mass) { next }
          if(!is.na(minPrecursorMass) & minPrecursorMass > mass) { next }
          if(is.na(maxSpectra) || maxSpectra > numberSpectraWritten) {
		parameterToCommand(settingsObject,fileName)
	 	numberSpectraWritten<-numberSpectraWritten+1
	  }
        }
        
        
      }
    }
    
  }
}
