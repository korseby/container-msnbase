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
if (length(args) < 4) {
    print("Error! No or not enough arguments given.")
    print("Usage: $0 dataMatrix sampleMetadata variableMetadata fileMSDialFormat")
    quit(save="no", status=1, runLast=FALSE)
}

# Data files
dataMatrix <- args[1]
sampleMetadata <- args[2]
variableMetadata <- args[3]
fileMSDialFormat <- args[4]



# ---------- convertToMSDialFormat ----------
dataMatrix$dataMatrix <- NULL

sampleMetadata <- t(sampleMetadata[,c(2,4,3)])

x <- ncol(variableMetadata)-2
y <- ncol(variableMetadata)-1
variableMetadata <- variableMetadata[,c(6,3,x,y)]

adductIonName <- !grepl(pattern="M", x=variableMetadata$adduct)
adductIon <- c("[M+H]+")
for(i in 1:nrow(variableMetadata)){
    if(adductIonName[i])
        variableMetadata[i,"adduct"] <- adductIon
}

# Start building Height table
heightTableMSDialFormat <- cbind(variableMetadata, dataMatrix)

# Delete isotopes
selIsotopes<-!grepl(pattern="M\\+", x=heightTableMSDialFormat$isotopes)
heightTableMSDialFormat <- heightTableMSDialFormat[selIsotopes,]

alignmentID <- c(1:nrow(heightTableMSDialFormat))
heightTableMSDialFormat <- cbind(alignmentID, heightTableMSDialFormat)

# Overwrite isotop-column with metabolite-column
metaboliteName <- rep("Unknown", nrow(heightTableMSDialFormat))
heightTableMSDialFormat[,4] <- metaboliteName

# Convert retention time from seconds to minutes
heightTableMSDialFormat[,"rt"] <- heightTableMSDialFormat[,"rt"]/60

# Names
sampleNames <- colnames(dataMatrix)
columnNames <- c("Alignment ID", "Average Rt(min)", "Average Mz", "Metabolite name", "Adduct ion name", sampleNames)
columnNames2 <- c("Class", "Type", "Injection order")

# Add first three rows
whiteSpaceMatrix <- matrix(data="", nrow=3, ncol=4)
firstThreeRows <- cbind(whiteSpaceMatrix, columnNames2, sampleMetadata)

# Build final Height table
colnames(heightTableMSDialFormat) <- columnNames
colnames(firstThreeRows) <- columnNames
rownames(heightTableMSDialFormat) <- NULL
rownames(firstThreeRows) <- NULL
heightTableMSDialFormat <- rbind(firstThreeRows, columnNames, heightTableMSDialFormat, stringsAsFactors=FALSE)
colnames(heightTableMSDialFormat) <- NULL

# Export file
dataSetName <- gsub(" ", "_", gsub(":", "_", Sys.time()))
fileMSDialFormat <- paste("Height_", dataSetName, sep="")
fileMSDialFormat <- paste(fileMSDialFormat, ".txt", sep="")

write.table(x=heightTableMSDialFormat, file=fileMSDialFormat, sep="\t", quote=FALSE, dec=".", row.names=FALSE)


