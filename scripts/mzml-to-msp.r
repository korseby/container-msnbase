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
    print("Usage: $0 working_dir plot_file polarity mzabs mzppm rtabs max.rt.range max.mz.range min.rt max.rt min.mz max.mz msms.intensity.threshold MergeSpectra")
    quit(save="no", status=1, runLast=FALSE)
}

# Working directory
setwd(dirname(args[1]))

# Data directory
mzml_dir <- dirname(args[1])

# Data file
mzml_files <- basename(args[1])

# Output plot file
plot_file <- args[2]

# MS1 variables
polarity <- as.character(args[3])
pol <- substr(x=polarity, start=1, stop=3)

# MS2 variables
mzabs <- as.numeric(args[4])                      # Absolute mass error (in seconds) used for merging MS/MS spectra
mzppm <- as.numeric(args[5])                      # ppm error used for merging MS/MS spectra
rtabs <- as.numeric(args[6])                      # Retention time error (in seconds) used for merging MS/MS spectra
max.rt.range <- as.numeric(args[7])               # Permitted retention time window (in seconds) of grouped MS1 precursors
max.mz.range <- as.numeric(args[8])               # Permitted m/z window of grouped MS1 precursors
min.rt <- as.numeric(args[9])                     # Minimum retention time for selected precursors
max.rt <- as.numeric(args[10])                    # Maximum retention time for selected precursors
min.mz <- as.numeric(args[11])                    # Minimum m/z value for selected precursors
max.mz <- as.numeric(args[12])                    # Maximum m/z value for selected precursors
msms.intensity.threshold <- as.numeric(args[13])  # Minimum intensity value for MS/MS peaks
MergeSpectra <- as.logical(args[14])              # Merge spectra enabled or disabled

# Preparations for plotting
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)



# ---------- Preparations ----------
# General variables
mzml_names <- NULL

# MS2 Matrices
msms_spectra <- NULL
msms_merged <- NULL



# ---------- Load mzML files ----------
f.load_mzml <- function() {
    # Load files
    #mzml_files <- list.files(mzml_dir, pattern="*.mzML", recursive=T, full.names=T)
    
    # Exclude blanks
    #mzml_files <- mzml_files[grep("MM8", mzml_files, invert=T)]
    #mzml_files <- mzml_files[grep("ACN", mzml_files, invert=T)]

    # Basenames of files without path and without extension
    #mzml_names <- gsub('(.*)\\..*', '\\1', gsub('( |-|,)', '.', basename(mzml_files)))
    mzml_names <- gsub('(.*)\\..*', '\\1', basename(mzml_files))
    
    # Return global variables
    mzml_files <<- mzml_files
    mzml_names <<- mzml_names
}



# ---------- Merge spectra group ----------
merge.spectra.group <- function(spectra, eppm, eabs, int.threshold) {
    if(is.na(spectra) || length(spectra) == 0) return(NULL)
    max.precursor.int <- -1
    rt.at.max.int <- NA
    start.rt <- attributes(spectra[[1]])$rt
    end.rt <- attributes(spectra[[1]])$rt
    mean.precursor.mz <- 0
    max.tic <- -1
    tandemms <- lapply(1:length(spectra), function(index) {
        # print(index)
        x<-spectra[[index]]
        # assign local variables
        cur.rt <- attributes(x)$rt
        cur.precursor.int <- attributes(x)$precursorIntensity
        cur.tic <- attributes(x)$tic
        cur.precursor.mz <- attributes(x)$precursorMz
        # assign global variables
        mean.precursor.mz <<- mean.precursor.mz + cur.precursor.mz
        max.tic <<- if(max.tic < cur.tic) cur.tic else max.tic
        if(max.precursor.int < cur.precursor.int) {
            max.precursor.int <<- cur.precursor.int
            rt.at.max.int <<- cur.rt
        }
        start.rt <<- if(start.rt > cur.rt) cur.rt else start.rt
        end.rt <<- if(end.rt < cur.rt) cur.rt else end.rt
        # start merging ms2 data
        a <- cbind(attributes(x)$mz, attributes(x)$intensity)
        colnames(a) <- c("mz", "intensity")
        return(a)
    })
    mean.precursor.mz <- mean.precursor.mz / length(spectra)
    
    # Filter intensity
    tandemms <- lapply(tandemms, function(spec) spec[spec[,"intensity"]>=int.threshold,])
    tandemrm <- sapply(1:length(tandemms), function(x) { if ( (dim(tandemms[[x]])[1] == 0) || is.null(dim(tandemms[[x]])[1]) ) x <- TRUE else x <- FALSE } )
    
    if ( (all(as.logical(tandemrm)) == TRUE) | (is.na(all(as.logical(tandemrm)))) )
        return(NULL)
    else
        tandemms <- tandemms[!tandemrm]
    
    # Process spectra
    peaks <- do.call(rbind, tandemms)
    g <- xcms:::mzClust_hclust(peaks[,"mz"], eppm=eppm, eabs=eabs)
    
    mz <- tapply (peaks[,"mz"], as.factor(g), mean)
    intensity <- tapply (peaks[,"intensity"], as.factor(g), max)
    
    if(length(mz) == 0)return(NULL);
    
    intensity <- intensity[order(mz)]
    mz <- mz[order(mz)]
    
    res <- .Call("Spectrum2_constructor_mz_sorted",
                 attributes(spectra[[1]])$msLevel, length(mz), rt.at.max.int, 
                 attributes(spectra[[1]])$acquisitionNum, 
                 attributes(spectra[[1]])$scanIndex, max.tic, mz,
                 intensity, attributes(spectra[[1]])$fromFile, attributes(spectra[[1]])$centroided, 
                 attributes(spectra[[1]])$smoothed, attributes(spectra[[1]])$polarity,
                 attributes(spectra[[1]])$merged, attributes(spectra[[1]])$precScanNum, 
                 mean.precursor.mz, max.precursor.int, attributes(spectra[[1]])$precursorCharge, attributes(spectra[[1]])$collisionEnergy,
                 TRUE, attributes(spectra[[1]])$.__classVersion__,
                 PACKAGE = "MSnbase")
    attributes(res)$start.rt <- start.rt
    attributes(res)$end.rt <- end.rt
    return(res)
}



# ---------- get.ppm.from.abs ----------
get.ppm.from.abs <- function(peak, ppm) {
    return ((peak / 1000000.0) * ppm)
}



# ---------- Collect spectra ----------
collect.spectra.lists <- function(spectra, mzabs, mzppm, rtabs) {
    a<-t(sapply(1:length(spectra), function(spectrum.id) {
        return(c(attributes(spectra[[spectrum.id]])$rt, attributes(spectra[[spectrum.id]])$precursorMz, spectrum.id))
    }))
    a <- a[order(a[,2], a[,1]),]
    index <- 1
    grouped.spectra.tmp <- list()
    while(index <= dim(a)[1]) {
        spectra.list <- list()
        spectra.list[[1]] <- spectra[[a[index,3]]]
        cur.abs <- get.ppm.from.abs(a[index, 2], mzppm * 0.5) + mzabs
        index <- index + 1
        while(index <= dim(a)[1]) {
            diff.mz <- abs(a[index - 1, 2] - a[index, 2])
            if(diff.mz <= cur.abs) {
                spectra.list[[length(spectra.list) + 1]] <- spectra[[a[index,3]]]
                index <- index + 1
            } else break
        }
        grouped.spectra.tmp[[length(grouped.spectra.tmp) + 1]] <- spectra.list
    }
    grouped.spectra <- list()
    sapply(1:length(grouped.spectra.tmp), function(spectrum.group.index) {
        spectrum.group <- grouped.spectra.tmp[[spectrum.group.index]]
        a<-t(sapply(1:length(spectrum.group), function(spectrum.index) {
            c(attributes(spectrum.group[[spectrum.index]])$rt, attributes(spectrum.group[[spectrum.index]])$precursorMz, spectrum.index)
        }))
        a <- matrix(a[order(a[,1], a[,2]),], ncol=3)
        index <- 1
        while(index <= dim(a)[1]) {
            spectra.list <- list()
            spectra.list[[1]] <- spectrum.group[[a[index,3]]]
            index <- index + 1
            while(index <= dim(a)[1]) {
                diff.rt <- abs(a[index - 1, 1] - a[index, 1])
                if(diff.rt <= rtabs) {
                    spectra.list[[length(spectra.list) + 1]] <- spectrum.group[[a[index,3]]]
                    index <- index + 1
                } else break
            }
            grouped.spectra[[length(grouped.spectra) + 1]] <<- spectra.list
        }
    })
    return(grouped.spectra)
}



# ---------- Filter grouped spectra ----------
filter.grouped.spectra <- function(grouped.spectra, max.rt.range, max.mz.range, min.rt, max.rt, min.mz, max.mz) {
    filtered.grouped.spectra <- list()
    sapply(grouped.spectra, function(grouped.spectrum) {
        mzs<-sapply(grouped.spectrum, function(x) attributes(x)$precursorMz) 
        rts<-sapply(grouped.spectrum, function(x) attributes(x)$rt)
        mean.mzs <- mean(mzs)
        mean.rts <- mean(rts)
        max.number.peaks <- max(sapply(grouped.spectrum, function(x) attributes(x)$peaksCount))
        if((max(mzs)-min(mzs)) <= max.mz.range & (max(rts)-min(rts)) <= max.rt.range & mean.rts < max.rt & mean.mzs < max.mz
           & mean.rts > min.rt & mean.mzs > min.mz & max.number.peaks > 0) {
            filtered.grouped.spectra[[length(filtered.grouped.spectra) + 1]] <<- grouped.spectrum
        }
    })
    return(filtered.grouped.spectra)
}



# ---------- Merge spectra ----------
merge.spectra <- function(spectra, mzabs, mzppm, rtabs, max.rt.range, max.mz.range, min.rt, max.rt, min.mz, max.mz, int.threshold = 0) {
    a <- collect.spectra.lists(spectra, mzabs, mzppm, rtabs)
    print(paste("Collected",length(a),"spectrum groups"))
    b <- filter.grouped.spectra(a, max.rt.range, max.mz.range, min.rt, max.rt, min.mz, max.mz)
    print(paste("Filtered",length(b),"spectrum groups"))
    merged.spectra <- sapply(1:length(b), function(x) {
        merge.spectra.group(b[[x]], mzppm * 10e-6, mzabs, int.threshold)
    })
    b<-Filter(Negate(is.null), merged.spectra)
    print(paste("Filtered",length(b),"spectrum groups with peak information"))
    return(b)
}



# ---------- Find MS/MS spectra ----------
f.ms2_find_spectra <- function() {
    # Find all MS/MS spectra in each sample
    msms_spectra <- list()
    msms_spectra <- readMSData(files=mzml_files, msLevel=2, verbose=TRUE)
    
    # Merge MS/MS spectra in each sample
    msms_merged <- list()
    if (MergeSpectra == TRUE) {
        msms_merged <- merge.spectra(msms_spectra, mzabs, mzppm, rtabs, max.rt.range, max.mz.range, min.rt, max.rt, min.mz, max.mz, msms.intensity.threshold)
    } else {
        msms_merged <- collect.spectra.lists(msms_spectra, mzabs, mzppm, rtabs)
        msms_merged <- lapply(msms_merged, FUN=function(x) { x <- x[[1]] })
    }

    # Return variables
    msms_spectra <<- msms_spectra
    msms_merged <<- msms_merged
}



# ---------- Plot MS/MS spectra ----------
f.ms2_plot_spectra <- function() {
    # For each profile
    msms_plot <- as.data.frame(matrix(unlist(lapply(msms_merged, FUN=function(x) { x <- data.frame(mz=x@precursorMz, rt=x@rt) } )), ncol=2))
    colnames(msms_plot) <- c("mz", "rt")
    msms_plot <- msms_plot[order(msms_plot$mz), ]
    
    pdf(plot_file, encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
    x_min <- floor(min(msms_plot[,"rt"]))
    x_max <- ceiling(max(msms_plot[,"rt"]))
    y_min <- floor(min(msms_plot[,"mz"]))
    y_max <- ceiling(max(msms_plot[,"mz"]))
    plot(x=msms_plot[,"rt"], y=msms_plot[,"mz"],
         xlim=c(x_min, x_max), ylim=c(y_min, y_max),
         main=mzml_names, xlab="retention time [rt]", ylab="merged m/z [mz]",
         pch=16, col=rainbow(n=nrow(msms_plot)))
    dev.off()
}



# ---------- Write out MS/MS spectra in MSP text format ----------
f.ms2_create_msp <- function() {
    msp_text <- NULL
    for (j in 1:length(msms_merged)) {
        NAME <- paste(mzml_names, ":", j, sep='')
        AlignmentID <- j
        RETENTIONTIME <- msms_merged[[j]]@rt
        PRECURSORMZ <- msms_merged[[j]]@precursorMz
        METABOLITENAME <- "Unknown"
        ADDUCTIONNAME <- "[M+H]+"
        NumPeaks <- msms_merged[[j]]@peaksCount
        Peaks <- NULL
        for (k in 1:length(msms_merged[[j]]@mz)) {
            Peaks <- c(Peaks, paste(msms_merged[[j]]@mz[k], msms_merged[[j]]@intensity[k], sep="\t") )
        }
        msp_text <- c(msp_text, paste("NAME:",NAME))
        msp_text <- c(msp_text, paste("AlignmentID:",AlignmentID))
        msp_text <- c(msp_text, paste("RETENTIONTIME:",RETENTIONTIME))
        msp_text <- c(msp_text, paste("PRECURSORMZ:",PRECURSORMZ))
        msp_text <- c(msp_text, paste("METABOLITENAME:",METABOLITENAME))
        msp_text <- c(msp_text, paste("ADDUCTIONNAME:",ADDUCTIONNAME))
        msp_text <- c(msp_text, paste("NumPeaks:",NumPeaks))
        msp_text <- c(msp_text, Peaks)
        msp_text <- c(msp_text, "")
    }
    cat(msp_text, file=paste(mzml_dir,"/",mzml_names,".msp",sep=""), sep="\n")
    print(paste("Creating ",mzml_dir,"/",mzml_names,".msp",sep=""))
}



# ---------- MAIN ----------
# Load files and define classes
f.load_mzml()
#f.sample_classes()

# MS2 Spectra Filtering and merging
f.ms2_find_spectra()
f.ms2_plot_spectra()
f.ms2_create_msp()

