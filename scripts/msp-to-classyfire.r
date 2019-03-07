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
if (length(args) < 2) {
    print("Error! No or not enough arguments given.")
    print("Usage: $0 mspfile classyfirefile")
    quit(save="no", status=1, runLast=FALSE)
}

# User variables
msp_file <- as.character(args[1])
classyfire_file <- as.numeric(args[2])



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



# ---------- ClassyFire ----------
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
cat(sep='\n', file=metfrag_tex_file, append=TRUE, c('\\begin{longtable}{ p{0.5cm} | p{1.5cm} | p{1.5cm} | p{3cm} | p{6cm} | p{1cm} }',
                                                    paste0("No", " & ", "RT", " & ", "MZ", " & ", "Mol. Formula", " & ", "Most abundand class", " & ", "Peaks", " \\\\"),
                                                    '\\hline'))
count <- 0
for (i in 1:nrow(metfrag_rtmz_unique)) {
    count <- count + 1
    metfrag_results_entries <- metfrag_results[c(which(metfrag_results$parentRT==metfrag_rtmz_unique[i,"parentRT"] & metfrag_results$parentMZ==metfrag_rtmz_unique[i,"parentMZ"])),]
    
    rt <- metfrag_results_entries$parentRT[1]
    mz <- metfrag_results_entries$parentMZ[1]
    molecular_formula <- metfrag_results_entries$MolecularFormula[1]
    most_abundant_class <- metfrag_results_entries$ClassyFireDirect
    most_abundant_class <- metfrag_results_entries$ClassyFireDirect[sort(table(most_abundant_class))[1]]
    mean_explained_peaks <- round(mean(metfrag_results_entries$NoExplPeaks, na.rm=TRUE), digits=0)
    number_peaks <- metfrag_results_entries$NumberPeaksUsed
    number_peaks <- metfrag_results_entries$NumberPeaksUsed[sort(table(number_peaks))[1]]
    
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
        cat(sep='\n', file=metfrag_tex_file, append=TRUE, c(paste0('			\\item[] \\textbf{Identifier:} ',metfrag_results_entries[j,"Identifier"]),
                                                            paste0('			\\item[] \\textbf{Synonyms:} ',metfrag_results_entries[j,"Synonyms"]),
                                                            paste0('			\\item[] \\textbf{Molecular Formula:} ',metfrag_results_entries[j,"MolecularFormula"]),
                                                            paste0('			\\item[] \\textbf{Primary compound class:} ', metfrag_results_entries[j,"ClassyFireDirect"]),
                                                            paste0('			\\item[] \\textbf{Alternative compound classes:} ', metfrag_results_entries[j,"ClassyFireAlternatives"]),
                                                            paste0('			\\item[] \\textbf{MetFusion Score:} ',metfrag_results_entries[j,"OfflineMetFusionScore"]),
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



