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
if (length(args) < 6) {
    print("Error! No or not enough arguments given.")
    print("Usage: $0 working_dir metfrag_results.csv metfrag_hits_limit metfrag_synonyms_method metfrag_classyfire_method metfrag_structure_method")
    quit(save="no", status=1, runLast=FALSE)
}

# Working directory
setwd(args[1])

# User variables
metfrag_results_file <- as.character(args[2])
metfrag_hits_limit <- as.numeric(args[3])
metfrag_synonyms_method <- as.logical(args[4]) # (True|False)
metfrag_classyfire_method <- as.logical(args[5]) # (True|False)
metfrag_structure_method <- as.character(args[6]) # (online|cdk)




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



# ---------- Fetch pictures of structures ----------
f.fetch_structures_pdf <- function(smiles) {
    #https://apps.ideaconsult.net/ambit2/depict/cdk/any?search=c1cc(c(c2c1C=C(CO2)c1ccc(cc1O)O)O)O&media=image/png
    #https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=C1%3DNC(%3DO)N%3DCN1&abbr=off&hdisp=bridgehead&showtitle=false&sma=cn&zoom=1.6&annotate=none
    metfrag_structure_url <- paste0('https://www.simolecule.com/cdkdepict/depict/bow/pdf?smi=', smiles, '&abbr=off&hdisp=bridgehead&showtitle=false&sma=cn&zoom=1.6&annotate=none')
    metfrag_structure_file <- paste0(tempfile(tmpdir="output"),'.pdf')
    result <- tryCatch(curl_fetch_disk(url=metfrag_structure_url, path=metfrag_structure_file), error=function(x) { x } )
    
    # Return pdf file location in which results are saved
    if (is.null(result$status_code)) {
        print("Error fetching structures file.")
        return(as.character(""))
    } else if ( (result$status_code >= 200) && (result$status_code <= 299) ) {
        return(as.character(metfrag_structure_file))
    } else {
        return(as.character(""))
    }
}



# ---------- Generate pictures of structures ----------
f.generate_structures_pdf <- function(smiles) {
    # Output PDF file
    pdf_file <- paste0(tempfile(tmpdir="output"),'.pdf')
    
    # Exec command
    suppressWarnings(
        cdk_output <- system2(command="java",
                              args=c('-jar','/usr/local/bin/cdk-smiles-to-pdf.jar',
                                     paste0("'",smiles,"'"),
                                     pdf_file),
                              stdout=TRUE, stderr=TRUE, wait=TRUE, timeout=60)
    )
    if (! is.null(attr(x=cdk_output, which="status"))) {
        print(paste0("Error: CDK exited with error code ", attr(x=cdk_output, which="status"), "."))
        return(as.character(""))
    }
    
    # Return pdf file location in which results are saved
    return(as.character(pdf_file))
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
metfrag_rtmz_unique <- metfrag_rtmz_unique[order(metfrag_rtmz_unique$parentRT), ]
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
metfrag_results[,"IUPACName"] <- gsub("\\$", "", metfrag_results[,"IUPACName"])
metfrag_results[,"IUPACName"] <- gsub("^", "\\\\( ", metfrag_results[,"IUPACName"])
metfrag_results[,"IUPACName"] <- gsub("$", " \\\\)", metfrag_results[,"IUPACName"])



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




# ---------- Structures as PDF ----------
metfrag_structures <- NULL

# Process each SMILES
print("Fetching structures...")
for (i in 1:nrow(metfrag_results)) {
    if (metfrag_structure_method == "online") {
        metfrag_structures <- c(metfrag_structures, f.fetch_structures_pdf(metfrag_results[i,"SMILES"]))
    } else { # "cdk"
        metfrag_structures <- c(metfrag_structures, f.generate_structures_pdf(metfrag_results[i,"SMILES"]))
    }
}

# Add column to metfrag_results
metfrag_results$StructurePDF <- gsub("output/", "", metfrag_structures)



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


# ---------- Generate TeX file ----------
# Process TeX file
metfrag_tex_file <- "output/metfrag_vis.tex"

# Write Header
cat(sep='\n', file=metfrag_tex_file, c('% Header',
'\\documentclass[9pt]{extarticle}',
'\\usepackage[T1]{fontenc}',
'\\usepackage[english]{babel}',
'\\usepackage{lmodern}',
'\\usepackage{tabularx}',
'\\usepackage{longtable}',
'\\usepackage[a4paper,total={170mm,257mm},left=20mm,top=20mm,margin=0.5in]{geometry}',
'\\usepackage{hyperref}',
'\\usepackage{graphicx}',
'\\begin{document}',
'\\setlength\\parindent{0pt}',
''))

# Write sample name
cat(sep='\n', file=metfrag_tex_file, append=TRUE, c(paste0('\\begin{center}\\section*{',metfrag_sample_name,'}\\end{center}'),
                                                    ''))

# Show summary at the beginning
cat(sep='\n', file=metfrag_tex_file, append=TRUE, c('\\begin{longtable}{ p{0.5cm} | p{1.5cm} | p{1.5cm} | p{3cm} | p{4.5cm} | p{3cm} }',
                                                    paste0("No", " & ", "RT", " & ", "MZ", " & ", "Mol. formula", " & ", "Most abundand class", " & ", "Explained peaks", " \\\\"),
                                                    '\\hline'))
count <- 0
for (i in 1:nrow(metfrag_rtmz_unique)) {
    count <- count + 1
    metfrag_results_entries <- metfrag_results[c(which(metfrag_results$parentRT==metfrag_rtmz_unique[i,"parentRT"] & metfrag_results$parentMZ==metfrag_rtmz_unique[i,"parentMZ"])),]
    
    rt <- metfrag_results_entries$parentRT[1]
    mz <- metfrag_results_entries$parentMZ[1]
    most_abundant_class <- metfrag_results_entries$ClassyFireDirect
    most_abundant_class <- metfrag_results_entries$ClassyFireDirect[sort(table(most_abundant_class))[1]]
    mean_explained_peaks <- round(mean(metfrag_results_entries$NoExplPeaks, na.rm=TRUE), digits=0)
    number_peaks <- metfrag_results_entries$NumberPeaksUsed
    number_peaks <- metfrag_results_entries$NumberPeaksUsed[sort(table(number_peaks))[1]]
    molecular_formula <- metfrag_results_entries$MolecularFormula[1]
    
    cat(sep='\n', file=metfrag_tex_file, append=TRUE, paste0(count, " & ", rt, " & ", mz, " & ", molecular_formula, " & ", most_abundant_class, " & ", paste0(mean_explained_peaks,"/",number_peaks), " \\\\"))
}

cat(sep='\n', file=metfrag_tex_file, append=TRUE, c('\\hline',
                                                    '\\end{longtable}',
                                                    '\\ \\\\',
                                                    '\\ \\\\'))

# Process MetFrag entries
print("Creating MetFrag TeX file...")
for (i in 1:nrow(metfrag_rtmz_unique)) {
    # Write entries for RT, MZ pairs
    metfrag_rtmz_unique[i,]
    metfrag_results_entries <- metfrag_results[c(which(metfrag_results$parentRT==metfrag_rtmz_unique[i,"parentRT"] & metfrag_results$parentMZ==metfrag_rtmz_unique[i,"parentMZ"])),]

    # Sort for RT, MZ entries
    cat(sep='\n', file=metfrag_tex_file, append=TRUE, c('\\noindent\\rule{\\textwidth}{0.4pt}',
                                                        paste0('\\subsection*{\\textbf{RT:} ',metfrag_results_entries[1,"parentRT"], ', \\textbf{MZ:} ', metfrag_results_entries[1,"parentMZ"], '}'),
                                                        '\\ \\\\',
                                                        '\\ \\\\',
                                                        ''))

    # Process each entry
    for (j in c(1:nrow(metfrag_results_entries))) {
        # Entry header
        cat(sep='\n', file=metfrag_tex_file, append=TRUE, c(paste0('% Entry: ',j),
                                                            '\\begin{minipage}{1\\textwidth}',
                                                            '    \\begin{minipage}{0.15\\textwidth}'))
        
        # Include structures
        if (nchar(metfrag_results_entries[j,"StructurePDF"]) > 2) {
            cat(sep='\n', file=metfrag_tex_file, append=TRUE, c(paste0('        \\includegraphics[width=4cm]{', metfrag_results_entries[j,"StructurePDF"], '}')))
        } else if (nchar(metfrag_results_entries[j,"StructurePDF"]) < 2) {
            cat(sep='\n', file=metfrag_tex_file, append=TRUE, c('        \\textbf{no structure available}'))
        } else if (file.info(metfrag_results_entries[j,"StructurePDF"])$size <= 2 ) {
            cat(sep='\n', file=metfrag_tex_file, append=TRUE, c('        \\textbf{no structure available}'))
        } else {
            cat(sep='\n', file=metfrag_tex_file, append=TRUE, c('        \\textbf{no structure available}'))
        }
        
        # Entry fill
        cat(sep='\n', file=metfrag_tex_file, append=TRUE, c('    \\end{minipage} \\hfill',
                                                            '    \\begin{minipage}{0.8\\textwidth}',
                                                            '        \\begin{itemize}'))
        
        # MetFrag items
        cat(sep='\n', file=metfrag_tex_file, append=TRUE, c(paste0('			\\item[] \\textbf{PubChem Identifier:} ',metfrag_results_entries[j,"Identifier"]),
                                                            paste0('			\\item[] \\textbf{Synonyms:} ',metfrag_results_entries[j,"Synonyms"]),
                                                            paste0('			\\item[] \\textbf{Molecular Formula:} ',metfrag_results_entries[j,"MolecularFormula"]),
                                                            paste0('			\\item[] \\textbf{Primary compound class:} ', metfrag_results_entries[j,"ClassyFireDirect"]),
                                                            paste0('			\\item[] \\textbf{Alternative compound classes:} ', metfrag_results_entries[j,"ClassyFireAlternatives"]),
                                                            paste0('			\\item[] \\textbf{Weighted MetFrag Score:} ',metfrag_results_entries[j,"Score"]),
                                                            paste0('			\\item[] \\textbf{Peaks explained:} ',metfrag_results_entries[j,"NoExplPeaks"], '/', metfrag_results_entries[j,"NumberPeaksUsed"]) ))
        
        # Entry footer
        cat(sep='\n', file=metfrag_tex_file, append=TRUE, c('        \\end{itemize}',
                                                            '    \\end{minipage}\\\\[0.4cm]',
                                                            '\\end{minipage}\\\\[0.8cm]',
                                                            ''))
    }
}

# TeX footer
cat(sep='\n', file=metfrag_tex_file, append=TRUE, c('% Footer',
                                                    '\\end{document}'))



# ---------- Generate PDF file ----------
#setwd(paste0(args[1],"/output"))
#suppressWarnings(
#    result <- system2(command="pdflatex", args="metfrag_vis.tex", stdout=TRUE, stderr=TRUE, #wait=TRUE, timeout=120)
#)
#if (! is.null(attr(x=result, which="status"))) {
#    cat(result, sep='\n')
#    print("Error generating PDF file.")
#    quit(status="400")
#}
#setwd(args[1])



