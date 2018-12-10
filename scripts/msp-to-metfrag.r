#!/usr/bin/env Rscript

# ---------- Preparations ----------
# Load libraries
# Load libraries
library(parallel)               # Detect number of cpu cores
library(foreach)                # For multicore parallel
library(doMC)                   # For multicore parallel
library(MSnbase)                # MS features
library(xcms)                   # Swiss army knife for metabolomics

# Setup R error handling to go to stderr
options(show.error.messages=F, error=function() { cat(geterrmessage(), file=stderr()); q("no",1,F) } )

# Set locales and encoding
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
loc <- Sys.setlocale(category="LC_ALL", locale="C")
options(encoding="UTF-8")

# Set options
options(stringAsfactors=FALSE, useFancyQuotes=FALSE)

# Multicore parallel
#nSlaves <- detectCores(all.tests=FALSE, logical=TRUE)
#registerDoMC(nSlaves)



# ---------- Arguments and user variables ----------
# Take in trailing command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 14) {
    print("Error! No or not enough arguments given.")
    print("Usage: $0 input.msp suspect.list output.txt minintthr minprop pretofrag fragtofrag polarity DatabaseSearchRelativeMassDeviation FragmentPeakMatchAbsoluteMassDeviation FragmentPeakMatchRelativeMassDeviation MetFragDatabaseType FilterExcludedElements FilterIncludedElements")
    quit(save="no", status=1, runLast=FALSE)
}

# Data file
msp_input_file <- args[1]
suspect_input_file <- args[2]
metfrag_output_file <- args[3]

# MSP variables
minintthr <- as.numeric(args[4])
minprop <- as.numeric(args[5])
pretofrag <- as.logical(args[6])
fragtofrag <- as.logical(args[7])

# MetFrag variables
polarity <- as.character(args[8])
pol <- substr(x=polarity, start=1, stop=3)
if (pol == "pos")
	IsPositiveIonMode <- "True"
else
	IsPositiveIonMode <- "False"

DatabaseSearchRelativeMassDeviation <- as.numeric(args[9])
FragmentPeakMatchAbsoluteMassDeviation <- as.numeric(args[10])
FragmentPeakMatchRelativeMassDeviation <- as.numeric(args[11])

MetFragDatabaseType <- as.character(args[12])

if (nchar(suspect_input_file) < 2) ScoreSuspectLists <- FALSE
	else if (suspect_input_file == "None") ScoreSuspectLists <- FALSE
	else ScoreSuspectLists <- TRUE

if (ScoreSuspectLists) {
	MetFragScoreTypes <- "FragmenterScore,OfflineMetFusionScore,SuspectListScore"
	MetFragScoreWeights <- "1.0,1.0,1.0"
} else {
	MetFragScoreTypes <- "FragmenterScore,OfflineMetFusionScore"
	MetFragScoreWeights <- "1.0,1.0"
}

MetFragPeakListReader <- "de.ipbhalle.metfraglib.peaklistreader.FilteredStringTandemMassPeakListReader"
MetFragCandidateWriter <- "CSV"
MaximumTreeDepth <- 2

FilterExcludedElements <- as.character(args[13])
FilterIncludedElements <- as.character(args[14])
MetFragPreProcessingCandidateFilter <- "UnconnectedCompoundFilter,IsotopeFilter,ElementInclusionOptionalFilter,ElementExclusionFilter"



# ---------- Preparations ----------
# General variables



# ---------- Parse MSP chunk ----------
parseMSP_chunk <- function(fileLines, minimumIntensityOfMaximalMS2peak, minimumProportionOfMS2peaks, neutralLossesPrecursorToFragments, neutralLossesFragmentsToFragments, offset = 0, progress = FALSE){
  
  ## LC-MS/MS entry:
  ## NAME: Unknown
  ## RETENTIONTIME: 3.215358
  ## PRECURSORMZ: 78.91963
  ## METABOLITENAME: 
  ## ADDUCTIONNAME: [M-H]-
  ## Num Peaks: 2
  ## 76.97093  754
  ## 76.98951  754
  ## 
  ## GC-MS additional properties:
  ## SCANNUMBER: 518
  ## MODELION: 59
  ## MODELIONHEIGHT: 924
  ## MODELIONAREA: 924
  ## INTEGRATEDHEIGHT: 924
  ## INTEGRATEDAREA: 924
  
  print("MS/MS file: Read file")
  
  #fileLines <- readLines(con = fileSpectra)
  numberOfFileLines <- length(fileLines)
  
  numberOfMS2PeaksOriginal <<- 0
  numberOfMS2PeaksWithNeutralLosses <<- 0
  numberOfTooHeavyFragments <<- 0
  numberOfMS2PeaksAboveThreshold <<- 0
  numberOfMS2PeaksBelowThreshold <<- 0
  numberOfSpectraDiscardedDueToNoPeaks      <<- 0
  numberOfSpectraDiscardedDueToMaxIntensity <<- 0
  numberOfSpectraDiscardedDueToTooHeavy <<- 0
  
  ## start with empty lines or not?
  endOfRecord <- TRUE
  if(numberOfFileLines > 0)
    if(nchar(trimws(fileLines[[1]])) > 0)
      endOfRecord <- FALSE
  
  ## check for pattern
  print("MS/MS file: Parse")
  isBI  	<- grepl(pattern = "^BEGIN IONS$",                    		x = fileLines)
  isName	<- grepl(pattern = "(^Name:)|(^NAME:)",               		x = fileLines)
  isNAme	<- grepl(pattern = "^NAME=",                           		x = fileLines)
  isNAME	<- grepl(pattern = "^TITLE=",                          		x = fileLines)
  isRT	  <- grepl(pattern = "^RETENTIONTIME:",						          x = fileLines)
  isRt	  <- grepl(pattern = "^retention time:",			  	          x = fileLines)
  isRts	  <- grepl(pattern = "^RTINSECONDS=",	     		  	          x = fileLines)
  isMZ	  <- grepl(pattern = "^PRECURSORMZ:",							          x = fileLines)
  isMz	  <- grepl(pattern = "^precursor m/z:",						          x = fileLines)
  isTotEm	<- grepl(pattern = "^total exact mass:",					        x = fileLines)
  isEMass	<- grepl(pattern = "(^exact mass:)|(^EXACT_MASS:)",       x = fileLines)
  isPMass	<- grepl(pattern = "^PEPMASS=",   							          x = fileLines)
  isE2Mass<- grepl(pattern = "^EXACTMASS=", 							          x = fileLines)
  isMetN	<- grepl(pattern = "^METABOLITENAME:",						        x = fileLines)
  isAddN	<- grepl(pattern = "(^ADDUCTIONNAME:)|(^Adductionname:)",	x = fileLines)
  isScanN	<- grepl(pattern = "^SCANNUMBER:",							          x = fileLines)
  isMIon	<- grepl(pattern = "^MODELION:",							            x = fileLines)
  isPrety	<- grepl(pattern = "^precursor type:",						        x = fileLines)
  isPretyp	<- grepl(pattern = "^PRECURSORTYPE:",						          x = fileLines)
  isNumP	<- grepl(pattern = "^Num Peaks:",							            x = fileLines)
  isPeak	<- grepl(pattern = "^\\d+(((\\.)|(,))\\d+)?[ \t]\\d+(((\\.)|(,))\\d+)?$",	x = fileLines)
  isCoCl	<- grepl(pattern = "^compound class:",						        x = fileLines)
  isInty	<- grepl(pattern = "^instrument type:",						        x = fileLines)
  isIntyp	<- grepl(pattern = "^INSTRUMENTTYPE:",						        x = fileLines)
  isIntype<- grepl(pattern = "^SOURCE_INSTRUMENT=",					        x = fileLines)
  isInt  	<- grepl(pattern = "^INSTRUMENT:",						            x = fileLines)
  isInchi	<- grepl(pattern = "(^InChI:)|(^InChI=)|(^INCHI=)|(^INCHI:)", x = fileLines)
  isInchiKey	<- grepl(pattern = "(^InChIKey:)|(^INCHIKEY:)|(^InChIKey=)|(^INCHIKEY=)|(^INCHIAUX=)",	      x = fileLines)
  isSmiles  	<- grepl(pattern = "(^SMILES:)|(^SMILES=)",		        x = fileLines)
  
  decimalDelimiterIsComma <- all(grepl(x = trimws(c(
    substring(text = fileLines[isMZ], first = nchar("PRECURSORMZ:") + 1), 
    substring(text = fileLines[isMz], first = nchar("precursor m/z:") + 1)
  )), pattern = "^(\\d+,\\d+$)|(^\\d+$)"))
  if(decimalDelimiterIsComma){
    fileLines_02 <- gsub(pattern = ",", replacement = ".", x = gsub(pattern = "\\.", replacement = "", x = fileLines))
  } else {
    fileLines_02 <- fileLines
  }
  
  ## extract
  suppressWarnings({
    parsedName	 <- 						      trimws(substring(text = fileLines, first = nchar("Name:") + 1))
    parsedNAme	 <- 						      trimws(substring(text = fileLines, first = nchar("NAME=") + 1))
    parsedNAME	 <- 						      trimws(substring(text = fileLines, first = nchar("TITLE=") + 1))
    parsedRT  	 <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("RETENTIONTIME:") + 1)))
    parsedRt     <- as.numeric(unlist(lapply(X = strsplit(x =   trimws(substring(text = fileLines_02, first = nchar("retention time:") + 1)), split = " "), FUN = function(x){if(length(x)==0) return(NA) else return(x[[1]])})))
    parsedRts  	 <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("RTINSECONDS=") + 1)))
    parsedMZ	   <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("PRECURSORMZ:") + 1)))
    parsedMz	   <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("precursor m/z:") + 1)))
    parsedTotEm  <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("total exact mass:") + 1)))
    parsedEMass  <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("exact mass:") + 1)))
    parsedE2Mass <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("EXACTMASS=") + 1)))
    parsedPMass  <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("PEPMASS=") + 1)))
    parsedMetN	 <- 						      trimws(substring(text = fileLines, first = nchar("METABOLITENAME:") + 1))
    parsedAddN	 <- 						      trimws(substring(text = fileLines, first = nchar("Adductionname:") + 1))
    parsedScanN  <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("SCANNUMBER:") + 1)))
    parsedMIon	 <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("MODELION:") + 1)))
    parsedPrety  <- 						      trimws(substring(text = fileLines, first = nchar("precursor type:") + 1))
    parsedPretyp <- 						      trimws(substring(text = fileLines, first = nchar("PRECURSORTYPE:") + 1))
    parsedNumP	 <- as.numeric(				trimws(substring(text = fileLines_02, first = nchar("Num Peaks:") + 1)))
    
    parsedTokensTmp <- strsplit(x = trimws(fileLines_02), split = "[ \t]")
    #parsedms2Peaks_mz <- as.numeric(unlist(lapply(X = parsedTokensTmp, FUN = function(x){head(x = x, n = 1)})))
    parsedms2Peaks_mz  <- as.numeric(unlist(lapply(X = parsedTokensTmp, FUN = function(x){if(length(x)<1) return(NA) else return(x[[1]])})))
    parsedms2Peaks_int <- as.numeric(unlist(lapply(X = parsedTokensTmp, FUN = function(x){if(length(x)<2) return(NA) else return(x[[2]])})))
    
    parsedCoCl	 <- 						      trimws(substring(text = fileLines, first = nchar("compound class:") + 1))
    parsedInty	 <- 						      trimws(substring(text = fileLines, first = nchar("instrument type:") + 1))
    parsedIntype <- 						      trimws(substring(text = fileLines, first = nchar("SOURCE_INSTRUMENT=") + 1))
    parsedIntyp  <- 						      trimws(substring(text = fileLines, first = nchar("INSTRUMENTTYPE:") + 1))
    parsedInt    <- 						      trimws(substring(text = fileLines, first = nchar("INSTRUMENT:") + 1))
    parsedInchi  <- 						      trimws(substring(text = fileLines, first = nchar("InChI:") + 1))
    parsedInchiKey  <- 						    trimws(substring(text = fileLines, first = nchar("InChIKey:") + 1))
    parsedSmiles    <- 						    trimws(substring(text = fileLines, first = nchar("SMILES:") + 1))
  })
  
  ## entry line intervals in file
  entryBorders   <- c(which(isName | isNAME | isBI), length(fileLines)+1)
  entryIntervals <- matrix(data = unlist(lapply(X = seq_len(length(entryBorders) - 1), FUN = function(x){c(entryBorders[[x]], entryBorders[[x+1]] - 1)})), nrow=2)
  
  numberOfSpectraOriginal <- length(entryBorders)
  
  ## do it
  print("MS/MS file: Assemble spectra")
  suppressWarnings(
    spectraList <- apply(X = entryIntervals, MARGIN = 2, FUN = function(x){
      #print(x)
      
      fileLines2      <- fileLines 	 		[x[[1]]:x[[2]]]
      fileLines_022   <- fileLines_02 	[x[[1]]:x[[2]]]
      parsedName2	 		<- parsedName	 		[x[[1]]:x[[2]]]
      parsedNAme2	 		<- parsedNAme	 		[x[[1]]:x[[2]]]
      parsedNAME2	 		<- parsedNAME	 		[x[[1]]:x[[2]]]
      parsedRT2  	 		<- parsedRT  	 		[x[[1]]:x[[2]]]
      parsedRt2     	<- parsedRt     	[x[[1]]:x[[2]]]
      parsedRts2     	<- parsedRts     	[x[[1]]:x[[2]]]
      parsedMZ2	 		  <- parsedMZ	 		  [x[[1]]:x[[2]]]
      parsedMz2	 		  <- parsedMz	 		  [x[[1]]:x[[2]]]
      parsedTotEm2  	<- parsedTotEm  	[x[[1]]:x[[2]]]
      parsedEMass2  	<- parsedEMass  	[x[[1]]:x[[2]]]
      parsedE2Mass2  	<- parsedE2Mass  	[x[[1]]:x[[2]]]
      parsedPMass2  	<- parsedPMass  	[x[[1]]:x[[2]]]
      parsedMetN2	 		<- parsedMetN	 		[x[[1]]:x[[2]]]
      parsedAddN2	 		<- parsedAddN	 		[x[[1]]:x[[2]]]
      parsedScanN2  	<- parsedScanN  	[x[[1]]:x[[2]]]
      parsedMIon2	 		<- parsedMIon	 		[x[[1]]:x[[2]]]
      parsedPrety2  	<- parsedPrety  	[x[[1]]:x[[2]]]
      parsedPretyp2  	<- parsedPretyp  	[x[[1]]:x[[2]]]
      parsedNumP2	 		<- parsedNumP	 		[x[[1]]:x[[2]]]
      parsedms2Peaks_mz2	<- parsedms2Peaks_mz	[x[[1]]:x[[2]]]
      parsedms2Peaks_int2	<- parsedms2Peaks_int	[x[[1]]:x[[2]]]
      parsedCoCl2	 		<- parsedCoCl	 		[x[[1]]:x[[2]]]
      parsedInty2	 		<- parsedInty	 		[x[[1]]:x[[2]]]
      parsedIntype2		<- parsedIntype		[x[[1]]:x[[2]]]
      parsedIntyp2		<- parsedIntyp		[x[[1]]:x[[2]]]
      parsedInt2	  	<- parsedInt  		[x[[1]]:x[[2]]]
      parsedInchi2  	<- parsedInchi  	[x[[1]]:x[[2]]]
      parsedInchiKey2 <- parsedInchiKey [x[[1]]:x[[2]]]
      parsedSmiles2  	<- parsedSmiles  	[x[[1]]:x[[2]]]
      
      isName2	 		<- isName	 		[x[[1]]:x[[2]]]
      isNAme2	 		<- isName	 		[x[[1]]:x[[2]]]
      isNAME2	 		<- isNAME	 		[x[[1]]:x[[2]]]
      isRT2  	 		<- isRT  	 		[x[[1]]:x[[2]]]
      isRt2     	<- isRt     	[x[[1]]:x[[2]]]
      isRts2     	<- isRts     	[x[[1]]:x[[2]]]
      isMZ2	 		  <- isMZ	 		  [x[[1]]:x[[2]]]
      isMz2	 		  <- isMz	 		  [x[[1]]:x[[2]]]
      isTotEm2  	<- isTotEm  	[x[[1]]:x[[2]]]
      isEMass2  	<- isEMass  	[x[[1]]:x[[2]]]
      isE2Mass2  	<- isE2Mass  	[x[[1]]:x[[2]]]
      isPMass2  	<- isPMass  	[x[[1]]:x[[2]]]
      isMetN2	 		<- isMetN	 		[x[[1]]:x[[2]]]
      isAddN2	 		<- isAddN	 		[x[[1]]:x[[2]]]
      isScanN2  	<- isScanN  	[x[[1]]:x[[2]]]
      isMIon2	 		<- isMIon	 		[x[[1]]:x[[2]]]
      isPrety2  	<- isPrety  	[x[[1]]:x[[2]]]
      isPretyp2  	<- isPretyp  	[x[[1]]:x[[2]]]
      isNumP2	 		<- isNumP	 		[x[[1]]:x[[2]]]
      isPeak2	    <- isPeak	    [x[[1]]:x[[2]]]
      isCoCl2	 		<- isCoCl	 		[x[[1]]:x[[2]]]
      isInty2	 		<- isInty	 		[x[[1]]:x[[2]]]
      isIntype2		<- isIntype		[x[[1]]:x[[2]]]
      isIntyp2		<- isIntyp		[x[[1]]:x[[2]]]
      isInt2  		<- isInt  		[x[[1]]:x[[2]]]
      isInchi2  	<- isInchi  	[x[[1]]:x[[2]]]
      isInchiKey2 <- isInchiKey [x[[1]]:x[[2]]]
      isSmiles2  	<- isSmiles  	[x[[1]]:x[[2]]]
      
      name <- NULL
      ms1Int <- NA
      rt <- NULL
      mz <- NULL
      metName <- "Unknown"
      adduct <- "Unknown"
      scanNumber <- NA
      quantMass <- NA
      peakNumber <- NA
      ms2Peaks_mz  <- vector(mode = "numeric")
      ms2Peaks_int <- vector(mode = "numeric")
      compoundClass <- "Unknown"
      instrumentType <- "Unknown"
      inchi <- ""
      inchiKey <- ""
      smiles <- ""
      
      if(any(isName2))
        ## name
        name        <- parsedName2	[[which(isName2)[[1]]]]
      if(any(isNAme2))
        ## name
        name        <- parsedNAme2	[[which(isNAme2)[[1]]]]
      if(any(isNAME2))
        ## name
        name        <- parsedNAME2	[[which(isNAME2)[[1]]]]
      if(any(isRT2))
        ## retention time
        rt          <- parsedRT2	[[which(isRT2)[[1]]]]
      if(any(isRt2) & any(is.null(rt), is.na(rt)))
        rt          <- parsedRt2	[[which(isRt2)[[1]]]]
      if(any(isRts2) & any(is.null(rt), is.na(rt)))
        rt          <- parsedRts2	[[which(isRts2)[[1]]]]
      if(any(isMZ2))
        ## precursor m/z
        mz          <- parsedMZ2	[[which(isMZ2)[[1]]]]
      if(any(isMz2)    & any(is.null(mz), is.na(mz), ifelse(test = is.null(mz), yes = FALSE, no = mz%%1==0)))
        mz          <- parsedMz2	[[which(isMz2)[[1]]]]
      if(any(isTotEm2) & any(is.null(mz), is.na(mz), ifelse(test = is.null(mz), yes = FALSE, no = mz%%1==0)) & any(isPrety2)){
        mzTmp <- parsedTotEm2[[which(isTotEm2)[[1]]]]
        pit   <- parsedPrety2[[which(isPrety2)[[1]]]]
        switch(pit,
               "[M-H]-" = { mz <- mzTmp -  1.008  },
               "[M-H]"  = { mz <- mzTmp -  1.008  },
               "[M+H]+" = { mz <- mzTmp +  1.008  },
               "[M+H]"  = { mz <- mzTmp +  1.008  },
               "[M+Na]+"= { mz <- mzTmp + 22.9898 },
               "[M+Na]" = { mz <- mzTmp + 22.9898 }#,
               #stop(paste("Unknown precursor ion type, pit, ")!", sep = ""))
        )
      }
      if(any(isEMass2) & any(is.null(mz), is.na(mz), ifelse(test = is.null(mz), yes = FALSE, no = mz%%1==0))){
      #if(any(isEMass2) & any(is.null(mz), is.na(mz), ifelse(test = is.null(mz), yes = FALSE, no = mz%%1==0)) & any(isPrety2)){
        mz <- parsedEMass2[[which(isEMass2)[[1]]]]
        #mzTmp <- parsedEMass2[[which(isEMass2)[[1]]]]
        #pit   <- parsedPrety2[[which(isPrety2)[[1]]]]
        #switch(pit,
        #       "[M-H]-" = { mz <- mzTmp -  1.008  },
        #       "[M-H]"  = { mz <- mzTmp -  1.008  },
        #       "[M+H]+" = { mz <- mzTmp +  1.008  },
        #       "[M+H]"  = { mz <- mzTmp +  1.008  },
        #       "[M+Na]+"= { mz <- mzTmp + 22.9898 },
        #       "[M+Na]" = { mz <- mzTmp + 22.9898 }#,
        #       #stop(paste("Unknown precursor ion type, pit, ")!", sep = ""))
        #)
      }
      if(any(isPMass2) & any(is.null(mz), is.na(mz), ifelse(test = is.null(mz), yes = FALSE, no = mz%%1==0))){
        if(!is.na(parsedPMass2[[which(isPMass2)[[1]]]])){
          mz  <- parsedPMass2[[which(isPMass2)[[1]]]]
        } else {
          tmp <- as.numeric(strsplit(x = trimws(substring(text = fileLines_022[[which(isPMass2)[[1]]]], first = nchar("PEPMASS=") + 1)), split = "\t")[[1]])
          mz  <- tmp[[1]]
          if(length(tmp) > 1)
            ms1Int <- tmp[[2]]
        }
      }
      if(any(isE2Mass2) & any(is.null(mz), is.na(mz), ifelse(test = is.null(mz), yes = FALSE, no = mz%%1==0))){
        mz <- parsedE2Mass2[[which(isE2Mass2)[[1]]]]
      }
      if(any(isMetN2))
        metName     <- parsedMetN2	[[which(isMetN2)[[1]]]]
      if(any(isAddN2)  & adduct == "Unknown")
        ## adduct
        adduct      <- parsedAddN2	[[which(isAddN2)[[1]]]]
      if(any(isPrety2) & adduct == "Unknown")
        adduct      <- parsedPrety2[[which(isPrety2)[[1]]]]
      if(any(isPretyp2) & adduct == "Unknown")
        adduct      <- parsedPretyp2[[which(isPretyp2)[[1]]]]
      if(any(isScanN2))
        scanNumber  <- parsedScanN2[[which(isScanN2)[[1]]]]
      if(any(isMIon2))
        quantMass  <- parsedMIon2	[[which(isMIon2)[[1]]]]
      if(any(isNumP2))
        ## #peaks
        peakNumber  <- parsedNumP2	[[which(isNumP2)[[1]]]]
      if(any(isPeak2)){
        ## MS2 peaks: "178.88669\t230"
        ms2Peaks_mz  <- parsedms2Peaks_mz2 [which(isPeak2)]
        ms2Peaks_int <- parsedms2Peaks_int2[which(isPeak2)]
      }
      if(any(isCoCl2))
        ## compound class
        compoundClass  <- parsedCoCl2	[[which(isCoCl2)[[1]]]]
      if(any(isInty2))
        ## instrument type
        instrumentType  <- parsedInty2	[[which(isInty2)[[1]]]]
      if(any(isIntype2))
        ## instrument type
        instrumentType  <- parsedIntype2	[[which(isIntype2)[[1]]]]
      if(any(isIntyp2))
        ## instrument type
        instrumentType  <- parsedIntyp2	[[which(isIntyp2)[[1]]]]
      if(all(any(isInt2), instrumentType %in% c("Unknown", "NA")))
        ## instrument
        instrumentType  <- parsedInt2	[[which(isInt2)[[1]]]]
      if(any(isInchi2))
        ## structure
        inchi  <- parsedInchi2 [[which(isInchi2)[[1]]]]
      if(any(isInchiKey2))
        ## structure
        inchiKey  <- parsedInchiKey2 [[which(isInchiKey2)[[1]]]]
      if(any(isSmiles2))
        ## structure
        smiles  <- parsedSmiles2 [[which(isSmiles2)[[1]]]]
      ## end of parsing
      
      if(is.null(rt))
        rt <- 0
      
      #if(is.null(mz))
      #  ## in case of gc
      #  mz <- max(ms2Peaks_mz)
      ms2Peaks_mz_original  <- ms2Peaks_mz
      ms2Peaks_int_original <- ms2Peaks_int
      
      numberOfMS2PeaksOriginal <<- numberOfMS2PeaksOriginal + length(ms2Peaks_mz)
      if(length(ms2Peaks_mz) == 0)
        numberOfSpectraDiscardedDueToNoPeaks <<- numberOfSpectraDiscardedDueToNoPeaks + 1
      
      ###################################################################
      ## filter fragments with mass greater than precursor
      numberOfTooHeavyFragmentsHere <- 0
      if(all(!is.null(mz), !is.na(mz))){
        tooHeavy <- ms2Peaks_mz > mz
        ms2Peaks_mz  <- ms2Peaks_mz [!tooHeavy]
        ms2Peaks_int <- ms2Peaks_int[!tooHeavy]
        numberOfTooHeavyFragmentsHere <- sum(tooHeavy)
        
        if(length(ms2Peaks_mz) == 0 & numberOfTooHeavyFragmentsHere > 0)
          numberOfSpectraDiscardedDueToTooHeavy <<- numberOfSpectraDiscardedDueToTooHeavy + 1
      }
      numberOfTooHeavyFragments <<- numberOfTooHeavyFragments + numberOfTooHeavyFragmentsHere
      
      ###################################################################
      ## filter for ms2 peak intensity
      peakNumber <- length(ms2Peaks_mz)
      if(peakNumber > 0){
        ###################################################################
        ## filter for ms2 peak intensity relative to maximum peak intensity in spectrum
        maximumIntensity <- max(ms2Peaks_int)
        if(maximumIntensity >= minimumIntensityOfMaximalMS2peak){
          ## spectrum is considered
          intensityThreshold <- maximumIntensity * minimumProportionOfMS2peaks
          fragmentsAboveThreshold <- ms2Peaks_int >= intensityThreshold
          
          ms2Peaks_mz  <- ms2Peaks_mz [fragmentsAboveThreshold]
          ms2Peaks_int <- ms2Peaks_int[fragmentsAboveThreshold]
          numberOfMS2PeaksAboveThreshold <<- numberOfMS2PeaksAboveThreshold + sum( fragmentsAboveThreshold)
          numberOfMS2PeaksBelowThreshold <<- numberOfMS2PeaksBelowThreshold + sum(!fragmentsAboveThreshold)
        } else {
          ## spectrum is not considered
          numberOfMS2PeaksBelowThreshold <<- numberOfMS2PeaksBelowThreshold + length(ms2Peaks_mz)
          numberOfSpectraDiscardedDueToMaxIntensity <<- numberOfSpectraDiscardedDueToMaxIntensity + 1
          ms2Peaks_mz  <- vector(mode = "numeric")
          ms2Peaks_int <- vector(mode = "numeric")
        }
      }
      
      ###################################################################
      ## normalize ms2 peaks to maximum = 1
      peakNumber <- length(ms2Peaks_mz)
      if(peakNumber > 0){
        max <- max(ms2Peaks_int)
        ms2Peaks_int <- ms2Peaks_int / max
      }
      
      ###################################################################
      ## add neutral losses
      if(peakNumber > 0){
        #################################
        ## neutral losses regarding the precursor
        if(all(!is.null(mz), !is.na(mz), neutralLossesPrecursorToFragments)){
          ms2PeaksNLPF_mz  <- ms2Peaks_mz - as.numeric(mz)
          ms2PeaksNLPF_int <- ms2Peaks_int
        } else {
          ms2PeaksNLPF_mz  <- vector(mode = "numeric")
          ms2PeaksNLPF_int <- vector(mode = "numeric")
        }
        #################################
        ## neutral losses amongst fragments
        if(neutralLossesFragmentsToFragments){
          m_mz  <- outer(X = ms2Peaks_mz,  Y = ms2Peaks_mz,  FUN = function(x,y){x-y})
          m_int <- outer(X = ms2Peaks_int, Y = ms2Peaks_int, FUN = function(x,y){(x+y) / 2})
          upper <- upper.tri(x = m_mz)
          ms2PeaksNLFF_mz  <- m_mz [upper]
          ms2PeaksNLFF_int <- m_int[upper]
        } else {
          ms2PeaksNLFF_mz  <- vector(mode = "numeric")
          ms2PeaksNLFF_int <- vector(mode = "numeric")
        }
        
        ms2Peaks_mz  <- c(ms2Peaks_mz,  ms2PeaksNLPF_mz,  ms2PeaksNLFF_mz)
        ms2Peaks_int <- c(ms2Peaks_int, ms2PeaksNLPF_int, ms2PeaksNLFF_int)
      }
      
      ###################################################################
      ## precursor mz
      #mz <- ifelse(test = !is.null(mz), yes = round(as.numeric(mz), digits = 4), no = ifelse(test = !is.null(scanNumber), yes = scanNumber, no = max(ms2Peaks_mz)))
      if(all(!is.null(mz), !is.na(mz))){
        mz <- round(as.numeric(mz), digits = 4)
      } else {
        if(!is.na(quantMass)){
          mz <- quantMass
        } else {
          if(!is.na(scanNumber)){
            mz <- scanNumber
          } else {
            mz <- max(ms2Peaks_mz)
          }
        }
      }
      
      ###################################################################
      ## string representation of spectrum
      spectrumString <- paste(ms2Peaks_mz_original, ms2Peaks_int_original, sep = " ", collapse = ";")
      
      ###################################################################
      ## built ms set
      spectrumItem <- list(
        name = name,
        ms1Int = ms1Int,
        #rt = round(as.numeric(rt), digits = 2),
        rt = rt,
        mz = mz,
        metName = metName,
        adduct = adduct,
        quantMass = quantMass,
        compoundClass = compoundClass,
        instrumentType = instrumentType,
        inchi = inchi,
        inchiKey = inchiKey,
        smiles = smiles,
        #peakNumber = as.numeric(peakNumber),
        peakNumber = length(ms2Peaks_mz),
        ms2Peaks_mz  = ms2Peaks_mz,
        ms2Peaks_int = ms2Peaks_int,
        spectrumString = spectrumString,
        entryInterval = x + offset
      )
      if(spectrumItem$peakNumber > 0){
        ## add
        numberOfMS2PeaksWithNeutralLosses <<- numberOfMS2PeaksWithNeutralLosses + spectrumItem$peakNumber
        return(spectrumItem)
      } else
        return(NULL)
    })
  )## suppressWarnings
  
  rm(
    isName,
    isNAME,
    isRT,
    isRt,
    isMZ,
    isMz,
    isTotEm,
    isEMass,
    isE2Mass,
    isPMass,
    isMetN,
    isAddN,
    isScanN,
    isMIon,
    isPrety,
    isNumP,
    isPeak,
    isCoCl,
    isInty,
    isInchi,
    isInchiKey,
    isSmiles,
    parsedName,
    parsedNAME,
    parsedRT,
    parsedRt,
    parsedMZ,
    parsedMz,
    parsedTotEm,
    parsedEMass,
    parsedE2Mass,
    parsedPMass,
    parsedMetN,
    parsedAddN,
    parsedScanN,
    parsedMIon,
    parsedPrety,
    parsedNumP,
    parsedTokensTmp,
    parsedms2Peaks_mz,
    parsedms2Peaks_int,
    parsedCoCl,
    parsedInty,
    parsedInchi,
    parsedInchiKey,
    parsedSmiles
  )
  
  print("MS/MS file: Box")
  
  ## remove NULL entries?
  spectraList[unlist(lapply(X = spectraList, FUN = is.null))] <- NULL
  
  numberOfSpectra <- length(spectraList)
  
  ## postprocess
  print(paste("MS/MS file postprocessing", sep = ""))
  #precursorMz <- vector(mode = "numeric", length = numberOfSpectra)
  #for(spectrumIdx in seq_len(length.out = numberOfSpectra))
  #  precursorMz[[spectrumIdx]] <- spectraList[[spectrumIdx]]$mz
  
  #precursorRt <- vector(mode = "numeric", length = numberOfSpectra)
  #for(spectrumIdx in seq_len(length.out = numberOfSpectra))
  #  precursorRt[[spectrumIdx]] <- spectraList[[spectrumIdx]]$rt
  
  precursorMz <- unlist(lapply(X = spectraList, FUN = function(x){
    if(is.null(x))  return(NULL)
    else            return(x$mz)
  }))
  precursorRt <- unlist(lapply(X = spectraList, FUN = function(x){
    if(is.null(x))  return(NULL)
    else            return(x$rt)
  }))
  
  print(paste("MS/MS file boxing", sep = ""))
  returnObj <- list()
  returnObj$fileSpectra <- NA
  returnObj$spectraList <- spectraList
  returnObj$numberOfSpectra <- numberOfSpectra
  returnObj$numberOfSpectraOriginal <- numberOfSpectraOriginal
  returnObj$numberOfMS2PeaksOriginal <- numberOfMS2PeaksOriginal
  returnObj$numberOfMS2PeaksWithNeutralLosses <- numberOfMS2PeaksWithNeutralLosses
  returnObj$numberOfMS2PeaksAboveThreshold <- numberOfMS2PeaksAboveThreshold
  returnObj$numberOfMS2PeaksBelowThreshold <- numberOfMS2PeaksBelowThreshold
  returnObj$numberOfTooHeavyFragments <- numberOfTooHeavyFragments
  returnObj$numberOfSpectraDiscardedDueToNoPeaks <- numberOfSpectraDiscardedDueToNoPeaks
  returnObj$numberOfSpectraDiscardedDueToMaxIntensity <- numberOfSpectraDiscardedDueToMaxIntensity
  returnObj$numberOfSpectraDiscardedDueToTooHeavy <- numberOfSpectraDiscardedDueToTooHeavy
  returnObj$precursorMz <- precursorMz
  returnObj$precursorRt <- precursorRt
  
  return(returnObj)
}



# ---------- Parse MSP ----------
parseMSP <- function(fileSpectra, minimumIntensityOfMaximalMS2peak, minimumProportionOfMS2peaks, neutralLossesPrecursorToFragments, neutralLossesFragmentsToFragments, progress = FALSE){
  fileLines <- readLines(con = fileSpectra)
  
  returnObj <- parseMSP_chunk(fileLines, minimumIntensityOfMaximalMS2peak, minimumProportionOfMS2peaks, neutralLossesPrecursorToFragments, neutralLossesFragmentsToFragments, progress)
  returnObj$fileSpectra <- fileSpectra
  
  return(returnObj)
}



# ---------- Convert MSP to MetFrag parameter ----------
f.ms2_msp_to_metfrag <- function() {
    # Read MSP
    parameterSet <- list()
    parameterSet$minimumIntensityOfMaximalMS2peak                  <- minintthr
    parameterSet$minimumProportionOfMS2peaks                       <- minprop
    parameterSet$neutralLossesPrecursorToFragments                 <- pretofrag
    parameterSet$neutralLossesFragmentsToFragments                 <- fragtofrag
    
    msp_file <- parseMSP(
        fileSpectra = msp_input_file, 
        minimumIntensityOfMaximalMS2peak = parameterSet$minimumIntensityOfMaximalMS2peak, 
        minimumProportionOfMS2peaks = parameterSet$minimumProportionOfMS2peaks, 
        neutralLossesPrecursorToFragments = parameterSet$neutralLossesPrecursorToFragments,
        neutralLossesFragmentsToFragments = parameterSet$neutralLossesFragmentsToFragments,
        progress = TRUE
    )
    
    # Write MetFrag parameter file
    metfrag_parameter_file <- NULL
    for (i in 1:msp_file$numberOfSpectra) {
        NeutralPrecursorMass <- msp_file$precursorMz[i]
        PrecursorIonType <- msp_file$spectraList[[i]]$adduct
        PeakListString <- gsub(" ", "_", msp_file$spectraList[[i]]$spectrumString)
        SampleName <- paste0("outputfolder/",
               i, "_",
               msp_file$spectraList[[i]]$rt, '_',
               msp_file$spectraList[[i]]$mz, "_",
               msp_file$spectraList[[i]]$ms1Int, "_",
               gsub(":.*", "", msp_file$spectraList[[i]]$name), ".mzML.txt")
        
        metfrag_parameter_line <- NULL
        metfrag_parameter_line <- paste0(metfrag_parameter_line,      "DatabaseSearchRelativeMassDeviation=", DatabaseSearchRelativeMassDeviation)
        metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "FragmentPeakMatchAbsoluteMassDeviation=", FragmentPeakMatchAbsoluteMassDeviation)
        metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "FragmentPeakMatchRelativeMassDeviation=", FragmentPeakMatchRelativeMassDeviation)
        metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "MetFragDatabaseType=", MetFragDatabaseType)
        if (ScoreSuspectLists) metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "ScoreSuspectLists=", suspect_input_file)
        metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "MetFragScoreTypes=", MetFragScoreTypes)
        metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "MetFragScoreWeights=", MetFragScoreWeights)
        metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "NeutralPrecursorMass=", NeutralPrecursorMass)
        metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "IsPositiveIonMode=", IsPositiveIonMode)
        metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "PrecursorIonType=", PrecursorIonType)
        metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "MetFragPeakListReader=", MetFragPeakListReader)
        metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "MetFragCandidateWriter=", MetFragCandidateWriter)
        metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "MaximumTreeDepth=", MaximumTreeDepth)
        metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "FilterExcludedElements=", FilterExcludedElements)
        metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "FilterIncludedElements=", FilterIncludedElements)
        metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "MetFragPreProcessingCandidateFilter=", MetFragPreProcessingCandidateFilter)
        metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "PeakListString=", PeakListString)
        metfrag_parameter_line <- paste0(metfrag_parameter_line, " ", "SampleName=", SampleName)
        
        metfrag_parameter_file <- c(metfrag_parameter_file, metfrag_parameter_line)
    }
    
    writeLines(text=metfrag_parameter_file, con=file(metfrag_output_file))
}



# ---------- MAIN ----------
# Read MSP spectra and create MetFrag parameter
f.ms2_msp_to_metfrag()



