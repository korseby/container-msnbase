#!/usr/bin/env Rscript

# ---------- Preparations ----------
# Load libraries
library(parallel)               # Detect number of cpu cores
library(foreach)                # For multicore parallel
library(doMC)                   # For multicore parallel
library(curl)                   # Downloading files with libcurl
library(jsonlite)               # For reading JSON files

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
if (length(args) < 3) {
    print("Error! No or not enough arguments given.")
    print("Usage: $0 metfrag_results.csv classyfire_results.csv classyfire_overview.pdf")
    quit(save="no", status=1, runLast=FALSE)
}

# User variables
metfrag_results_file <- as.character(args[1])
metfrag_output_file <- as.character(args[2])
metfrag_plot_file <- as.character(args[3])
metfrag_hits_limit <- as.numeric(1)
metfrag_synonyms_method <- as.logical(TRUE)
metfrag_classyfire_method <- as.logical(TRUE)



# ---------- Fetch synonyms for a SMILES ----------
f.synonyms_smiles <- function(smiles) {
    # Check SMILES
    if ( (is.na(smiles)) || (smiles == "") ) {
        return(as.character(""))
    }
    
    # Fetch synonyms
    metfrag_synonyms_url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/", smiles, "/synonyms/JSON")
    result <- tryCatch(fromJSON(metfrag_synonyms_url), error=function(x) { x } )
    
    # No synonyms found
    if (is.null(result$InformationList)) {
        return(as.character(""))
    }
    
    # Synonyms
    metfrag_synonyms <- gsub("[^[:lower:][:upper:][:digit:]/\\-]", "", result$InformationList$Information$Synonym[[1]])
    if (length(metfrag_synonyms) > 10)
        metfrag_synonyms <- metfrag_synonyms[c(1:10)]
    return(paste(metfrag_synonyms, sep="", collapse='; '))
}



# ---------- Classyfire InChI ----------
f.classyfire_inchi <- function(inchi) {
    # Example: http://classyfire.wishartlab.com/queries/3339010
    inchi <- as.character(inchi)
    
    if (is.na(inchi)) {
        print(paste0("        Error: InChI empty. Skipping ClassyFire."))
        return(list(NULL))
    } else if (nchar(inchi) < 2) {
        print(paste0("        Error: InChI empty. Skipping ClassyFire."))
        return(list(NULL))
    } else {
        # Send request to ClassyFire
        suppressWarnings(
            classyfire_output <- system2(command="curl",
                                         args=c('--header "Content-Type: application/json"',
                                                '--include --silent --request POST',
                                                "--data '{\"label\":\"metfrag\",\"query_input\":\"",inchi,"\",\"query_type\":\"STRUCTURE\"}'",
                                                '"http://classyfire.wishartlab.com/queries.json"'),
                                         stdout=TRUE, stderr=TRUE, wait=TRUE, timeout=60)
        )
        if (! is.null(attr(x=classyfire_output, which="status"))) {
            print(paste0("        Error: ClassyFire exited with error code ", attr(x=classyfire_output, which="status"), "."))
            return(list(NULL))
        }
        
        if ( ! ( (grepl("HTTP/1.1 201 Created", classyfire_output)) || (! grepl("HTTP/1.1 4", classyfire_output)) || (! grepl("HTTP/1.1 5", classyfire_output)) ) ) {
            print("        Error: ClassyFire exited with error.")
            return(list(NULL))
        }
        
        # Query ID
        classyfire_id <- strsplit(x=strsplit(x=classyfire_output[length(classyfire_output)],split=",")[[1]][[1]],split=":")[[1]][[2]]
        
        # Wait for results to be complete
        print(paste("        Request-URL: http://classyfire.wishartlab.com/queries/",classyfire_id,sep=''))
        wait <- 1
        while (TRUE) {
            # Fetch ClassyFire JSONs
            classyfire_results <- tryCatch(jsonlite::fromJSON(txt=url(paste('http://classyfire.wishartlab.com/queries/',classyfire_id,'.json',sep=''))), error=function(x) { NA } )
            
            if (length(classyfire_results) <= 1) {
                break
            } else if ( (classyfire_results$classification_status == "Done") || 
                        (classyfire_results$classification_status == "In progress") ||
                        (wait >= 6) ) {
                break
            }
            print(paste0("        Status: ", classyfire_results$classification_status, ". Polling ClassyFire for results to be complete..."))
            wait <- wait + 1
            Sys.sleep(time=10)
        }
        if (length(classyfire_results) <= 1) {
            print("        Error fetching JSON file.")
            classyfire_results <- NULL
        } else if ( (grepl("Invalid Request", classyfire_results$classification_status)) ||
                    (grepl("Internal Server Error", classyfire_results$classification_status)) ||
                    (grepl("In Queue", classyfire_results$classification_status)) ) {
            print(paste0("        Error: ", classyfire_results$classification_status, "!!!"))
            classyfire_results <- NULL
        }
        
        # Return results
        return(list(classyfire_results$entities))
    }
}



# ---------- Process MetFrag Table ----------
# Read MetFrag results file
metfrag_results_raw <- read.table(file=metfrag_results_file, quote='\"', sep=',', header=TRUE, stringsAsFactors=TRUE, fill=TRUE)

# Only take the entries with the highest scores
metfrag_results <- metfrag_results_raw[0, ]
metfrag_rtmz_unique <- unique(metfrag_results_raw[,c("parentRT","parentMZ")])
metfrag_rtmz_unique <- metfrag_rtmz_unique[order(metfrag_rtmz_unique$parentMZ), ]
for (i in 1:nrow(metfrag_rtmz_unique)) {
    temp <- metfrag_results_raw[c(which(metfrag_results_raw$parentRT==metfrag_rtmz_unique[i,"parentRT"] & metfrag_results_raw$parentMZ==metfrag_rtmz_unique[i,"parentMZ"])),]
    temp <- temp[order(temp$Score, decreasing=TRUE),]
    temp <- temp[c(1:metfrag_hits_limit),]
    metfrag_results <- rbind(metfrag_results, temp)
}



# ---------- Normalize names ----------
# Normalize MetFrag sample name
metfrag_sample_name <- metfrag_results_raw[1,"fileName"]
metfrag_sample_name <- gsub("\\..*", "", metfrag_sample_name, ignore.case=TRUE)
metfrag_sample_name <- gsub("[^0-9A-Za-z///' ]", " ", metfrag_sample_name, ignore.case=TRUE)

# Normalize IUPACName for TeX
#metfrag_results[,"IUPACName"] <- gsub("\\$", "", metfrag_results[,"IUPACName"])
#metfrag_results[,"IUPACName"] <- gsub("^", "\\\\( ", metfrag_results[,"IUPACName"])
#metfrag_results[,"IUPACName"] <- gsub("$", " \\\\)", metfrag_results[,"IUPACName"])



# ---------- Fetch synonyms ----------
if (metfrag_synonyms_method == TRUE) {
    metfrag_synonyms <- NULL
    
    # Process each InChI
    print("Fetching synonyms...")
    for (i in 1:nrow(metfrag_results)) {
        metfrag_synonyms <- c(metfrag_synonyms, f.synonyms_smiles(metfrag_results[i,"SMILES"]))
    }
    
    # Add column to metfrag_results
    metfrag_results$Synonyms <- metfrag_synonyms
} else {
    metfrag_results$Synonyms <- ""
}




# ---------- ClassyFire ----------
if (metfrag_classyfire_method == TRUE) {
    # Create empty list
    metfrag_classyfire <- list()
    
    # Process each InChI
    for (i in 1:nrow(metfrag_results)) {
        print(paste("Classifying MetFrag result",i,"of",nrow(metfrag_results),"..."))
        
        metfrag_classyfire <- append(metfrag_classyfire, as.list(f.classyfire_inchi(metfrag_results[i,"InChI"])))
    }
    
    # Add classyfire columns to metfrag_results
    metfrag_results$ClassyFireDirect <- as.character(lapply(X=metfrag_classyfire, FUN = function(x) { x <- x$direct_parent$name } ))
    metfrag_results$ClassyFireAlternatives <- as.character(lapply(X=metfrag_classyfire, FUN = function(x) { paste(as.data.frame(x$alternative_parents)$name, sep="", collapse="; ") } ))
} else {
    metfrag_results$ClassyFireDirect <- ""
    metfrag_results$ClassyFireAlternatives <- ""
}



# ---------- Save results as CSV ----------
# Write MetFrag results file
write.csv(x=metfrag_results, file=metfrag_output_file, quote=TRUE, row.names=FALSE, na="\"\"")



# ---------- Statistics on classes ----------
# Check whether there are known compounds
print(paste(length(which(metfrag_results$Synonyms != "")), "identified compounds."))

# Determine most abundant primary classes
primary_classes <- data.frame(classes=as.character(unique(metfrag_results[(metfrag_results$ClassyFireDirect != "NULL"), "ClassyFireDirect"])), frequency=0)
for (i in 1:length(primary_classes$classes)) primary_classes[i,"frequency"] <- length(which(metfrag_results$ClassyFireDirect == primary_classes[i,"classes"]))
primary_classes <- primary_classes[order(primary_classes$frequency, decreasing=TRUE),]

# Determine most abundant subordinate classes
temp <- data.frame(classes=unlist(strsplit(x=paste(metfrag_results$ClassyFireAlternatives, collapse="; "), split="; ")), frequency=0)
if ( length(which(temp$classes == "")) > 0) temp <- temp[-which(temp$classes == ""),]
secondary_classes <- data.frame(classes=unique(temp$classes), frequency=0)
for (i in 1:length(secondary_classes$classes)) secondary_classes[i,"frequency"] <- length(which(temp$classes == secondary_classes[i,"classes"]))
secondary_classes <- secondary_classes[order(secondary_classes$frequency, decreasing=TRUE),]

pdf(file=metfrag_plot_file, encoding="ISOLatin1", pointsize=10, width=20, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(14,4,4,1), oma=c(0,0,0,0), cex.axis=0.8, cex=0.9)
barplot(primary_classes[1:20,"frequency"], names.arg=primary_classes[1:20,"classes"], las=3, main="Most abundand primary classes")
barplot(secondary_classes[1:40,"frequency"], names.arg=secondary_classes[1:40,"classes"], las=3, main="Most abundand subordinate classes")
dev.off()



