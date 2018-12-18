#!/usr/bin/env Rscript

# ---------- Preparations ----------
# Load libraries
# Load libraries
library(parallel)               # Detect number of cpu cores
library(foreach)                # For multicore parallel
library(doMC)                   # For multicore parallel

# Setup R error handling to go to stderr
options(show.error.messages=F, error=function() { cat(geterrmessage(), file=stderr()); q("no",1,F) } )

# Set locales and encoding
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
loc <- Sys.setlocale(category="LC_ALL", locale="C")
options(encoding="UTF-8")

# Set options
options(stringAsfactors=FALSE, useFancyQuotes=FALSE)

# Multicore parallel
nSlaves <- detectCores(all.tests=FALSE, logical=TRUE)
registerDoMC(nSlaves)



# ---------- MAIN ----------
print(paste("Number of detected CPU cores:",as.character(nSlaves)))

