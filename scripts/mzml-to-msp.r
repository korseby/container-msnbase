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
#options(show.error.messages=F, error=function() { cat(geterrmessage(), file=stderr()); q("no",1,F) } )

# Set locales and encoding
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
loc <- Sys.setlocale(category="LC_ALL", locale="C")
options(encoding="UTF-8")

# Set options
options(stringAsfactors=FALSE, useFancyQuotes=FALSE)

# Multicore parallel
nSlaves <- detectCores(all.tests=FALSE, logical=TRUE)
registerDoMC(nSlaves)



# ---------- Arguments and user variables ----------
# Take in trailing command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
    print("Error! No or not enough arguments given.")
    print("Usage: $0 working_dir")
    quit(save="no", status=1, runLast=FALSE)
}

# Working directory
setwd(args[1])

# Data directory
mzml_dir <- args[1]

# MS1 variables
polarity <- "positive"
pol <- substr(x=polarity, start=1, stop=3)
ppm <- 35
peakwidth <- c(4,21)
snthresh <- 10
prefilter <- c(5,50)
fitgauss <- FALSE
verbose.columns <- FALSE
mzwidth <- 0.01
minfrac <- 0.5
bwindow <- 4

# MS2 variables
mzabs <- 0.01                         # Absolute mass error (in seconds) used for merging MS/MS spectra
mzppm <- 5                            # ppm error used for merging MS/MS spectra
rtabs <- 5                            # Retention time error (in seconds) used for merging MS/MS spectra
max.rt.range <- 20                    # Permitted retention time window (in seconds) of grouped MS1 precursors
max.mz.range <- 0.01                  # Permitted m/z window of grouped MS1 precursors
min.rt <- 10                          # Minimum retention time for selected precursors
max.rt <- 1020                        # Maximum retention time for selected precursors
min.mz <- 50                          # Minimum m/z value for selected precursors
max.mz <- 1500                        # Maximum m/z value for selected precursors
msms.intensity.threshold <- 100       # Minimum intensity value for MS/MS peaks

# Preparations for plotting
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)



# ---------- Preparations ----------
# General variables
mzml_files <- NULL
mzml_names <- NULL

# Specific variables
species <- NULL
seasons <- NULL
seasonal_species <- NULL
spesearep <- NULL
species_names <- NULL
species_colors <- NULL
species_symbols <- NULL
seasons_names <- NULL
seasons_colors <- NULL
seasons_symbols <- NULL
species_samples_colors <- NULL
seasons_samples_colors <- NULL
species_samples_symbols <- NULL
seasons_samples_symbols <- NULL
samples <- NULL
presence <- 12 * (10/12)

# MS1 Matrices
peak_xset <- NULL
peak_xcam <- NULL
diff_list <- NULL
peak_list <- NULL
feat_list <- NULL
bina_list <- NULL
uniq_list <- NULL

# MS2 Matrices
msms_spectra <- NULL
msms_merged <- NULL



# ---------- Load mzML files ----------
f.load_mzml <- function() {
    # Load files
    mzml_files <- list.files(mzml_dir, pattern="*.mzML", recursive=T, full.names=T)
    
    # Exclude blanks
    #mzml_files <- mzml_files[grep("MM8", mzml_files, invert=T)]
    #mzml_files <- mzml_files[grep("ACN", mzml_files, invert=T)]

    # Basenames of files without path and without extension
    mzml_names <- gsub('(.*)\\..*', '\\1', gsub('( |-|,)', '.', basename(mzml_files)))
    
    # Return global variables
    mzml_files <<- mzml_files
    mzml_names <<- mzml_names
}



# ---------- Define sample classes ----------
f.sample_classes <- function() {
    # Sample classes: species
    species <<- as.factor(sapply(strsplit(as.character(mzml_names), "_"), function(x) {
        nam <- x[4];
        nam;
    }))
    
    # Sample classes: seasons
    seasons <<- as.factor(sapply(strsplit(as.character(mzml_names), "_"), function(x) {
        nam <- x[3];
        nam;
    }))
    
    # Sample classes: seasonal species
    seasonal_species <<- as.factor(sapply(strsplit(as.character(mzml_names), "_"), function(x) {
        se <- x[3];
        sp <- x[4];
        nam <- paste(se, '_', sp, sep='')
        nam;
    }))
    
    # Sample classes: unique species-seasons-replicate
    spesearep <<- as.factor(sapply(strsplit(as.character(mzml_files), "_"), function(x) {
        se <- x[3];
        sp <- x[4];
        nam <- paste(sp, '_', se, sep='')
        nam;
    }))
    spesearep <<- make.names(as.character(spesearep), unique=TRUE)
    
    # Define species names, colors, symbols
    species_names <<- levels(species)
    species_colors <<- c("yellowgreen", "mediumseagreen", "darkorange1", "firebrick3", "darkolivegreen4", "dodgerblue4", "chocolate", "darkviolet", "darkkhaki")
    species_symbols <<- c(15, 16, 0, 1, 17, 8, 2, 5, 18)
    
    # Define seasons names, colors, symbols
    seasons_names <<- c("summer", "autumn", "winter", "spring")
    seasons_colors <<- c("darkgoldenrod3", "firebrick3", "deepskyblue3", "chartreuse3")
    seasons_symbols <<- c(1, 15, 16, 0)
    
    # Define samples colors
    species_samples_colors <<- sapply(species, function(x) { x <- species_colors[which(x==species_names)] } )
    seasons_samples_colors <<- sapply(seasons, function(x) { x <- seasons_colors[which(x==seasons_names)] } )
    
    # Define samples symbols
    species_samples_symbols <<- sapply(species, function(x) { x <- species_symbols[which(x==species_names)] } )
    seasons_samples_symbols <<- sapply(seasons, function(x) { x <- seasons_symbols[which(x==seasons_names)] } )
    
    # Define number of samples
    samples <<- length(species)/length(species_names)
}



# ---------- Find MS/MS spectra ----------
f.ms2_find_spectra <- function() {
    # Find all MS/MS spectra in each sample
    msms_spectra <- list()
    msms_spectra <- foreach(i=1:length(mzml_files)) %dopar% {
        readMSData(files=mzml_files[i], msLevel=2, verbose=TRUE)
    }
    
    # Merge MS/MS spectra in each sample
    msms_merged <- list()
    msms_merged <- foreach(i=1:length(mzml_files)) %dopar% {
        #print(paste("Merging MS/MS spectras in",msms_spectra[[i]]@phenoData@data$sampleNames,"..."))
        merge.spectra(msms_spectra[[i]], mzabs, mzppm, rtabs, max.rt.range, max.mz.range, min.rt, max.rt, min.mz, max.mz, msms.intensity.threshold)
    }
    
    
    # Plot number of MS/MS spectra in each sample
    #msms_num_spectra <- data.frame(species=species, sample=mzml_names, nspectra=0)
    #for (i in 1:length(msms_merged))
    #    msms_num_spectra[i,"nspectra"] <- length(msms_merged[[i]])
    
    # For species
    #boxplot(nspectra ~ species, data=msms_num_spectra, col=species_colors, names=NA,
    #        main="Number of acquired Auto-MS spectra per species", xlab="Species", ylab="Number of spectra")
    #text(1:length(species_names), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=species_names, xpd=TRUE, cex=0.9)
    
    # For seasons
    # R-bug: prevent sorting when using formula
    #model_boxplot <- data.frame(summer=msms_num_spectra$nspectra[seasons=="summer"],
    #                            autumn=msms_num_spectra$nspectra[seasons=="autumn"],
    #                            winter=msms_num_spectra$nspectra[seasons=="winter"],
    #                            spring=msms_num_spectra$nspectra[seasons=="spring"])
    #boxplot(model_boxplot, data=msms_num_spectra, col=seasons_colors,
    #        main="Number of acquired Auto-MS spectra in the seasons", xlab="Seasons", ylab="Number of spectra")
    
    # Return variables
    msms_spectra <<- msms_spectra
    msms_merged <<- msms_merged
}



# ---------- Write out MS/MS spectra in MSP text format ----------
f.ms2_create_msp <- function() {
    for (i in 1:length(msms_merged)) {
        msp_text <- NULL
        for (j in 1:length(msms_merged[[i]])) {
            NAME <- paste(mzml_names[i], ":", j, sep='')
            AlignmentID <- j
            RETENTIONTIME <- msms_merged[[i]][[j]]@rt
            PRECURSORMZ <- msms_merged[[i]][[j]]@precursorMz
            METABOLITENAME <- "Unknown"
            ADDUCTIONNAME <- "[M+H]+"
            NumPeaks <- msms_merged[[i]][[j]]@peaksCount
            Peaks <- NULL
            for (k in 1:length(msms_merged[[i]][[j]]@mz)) {
                Peaks <- c(Peaks, paste(msms_merged[[i]][[j]]@mz[k], msms_merged[[i]][[j]]@intensity[k], sep="\t") )
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
        cat(msp_text, file=paste(mzml_dir,"/../msp/",mzml_names[i],".msp",sep=""), sep="\n")
    }
}



# ---------- MAIN ----------
# Load files and define classes
f.load_mzml()
#f.sample_classes()

# MS2 Spectra Filtering and merging
f.ms2_find_spectra()
f.ms2_create_msp()


