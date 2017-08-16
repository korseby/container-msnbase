# MSnbase
Version: 2.0.2

## Short Description
Basic plotting, data manipulation and processing of MS-based Proteomics data.

## Description

MSnbase aims are providing a reproducible research framework to proteomics data analysis. It should allow researcher to easily mine mass spectrometry data, explore the data and its statistical properties and visually display these.

MSnbase also aims at being compatible with the infrastructure implemented in Bioconductor, in particular Biobase. As such, classes developed specifically for proteomics mass spectrometry data are based on the eSet and ExpressionSet classes. The main goal is to assure seamless compatibility with existing meta data structure, accessor methods and normalisation techniques.

## Key features

- Processing LC/MS and GC/MS data

## Approaches

- Proteomics / Untargeted
- Proteomics / Targeted
- Metabolomics / Untargeted
- Metabolomics / Targeted

## Instrument Data Types

- MS / LC-MS
- MS / GC-MS

## Tool Authors

Laurent Gatto
Guangchuang Yu
Samuel Wieczorek
Vasile-Cosmin Lazar
Vladislav Petyuk
Thomas Naake
Richie Cotton
Martina Fisher
Johannes Rainer
Sebastian Gibb

## Container Contributors

- [Kristian Peters](https://github.com/korseby) (IPB-Halle)

## Git Repository

- https://github.com/lgatto/MSnbase

## Publications

- Gatto L, Lilley K (2012): MSnbase - an R/Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation. Bioinformatics 28: 288-289.

## Installation 

```bash
docker build -t msnbase .
```
Alternatively, pull from repo:
```bash
docker pull container-registry.phenomenal-h2020.eu/phnmnl/msnbase
```

## Usage Instructions
On a PhenoMeNal Cloud Research Environment Galaxy environment, go to MS tool category or type msnbase-read-msms in the search tools text field, and then click on msnbase-read-msms and select a mzML file containing MS/MS information from the history, then press run.

