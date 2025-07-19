# SMTrackR

## Visualizing protein-DNA binding on individual sequenced DNA molecules.

Single-molecule assays like NOMe-seq, and dSMF are superior to DNase-seq and
ATAC-seq as they do not nibble down DNA. Thus, they enable quantification of
all three i.e., protein-free, Transcription Factor-bound and
histone-complex-bound states. But a user-friendly tool to visualize and
quantify such states is lacking. Here, we present, SMTrackR, a Bioconductor
package to visualize protein-DNA binding states on individual sequenced DNA
molecules. SMTrackR queries the single-molecule footprint database we built and
hosted at Galaxy Server. It comprises of BigBed files generated from NOMe-seq
and dSMF datasets. SMTrackR exploits UCSC REST API to query a BigBed file and
plot footprint heatmap categorized in different binding states as well as
report their occupancies.

List of tracks are available [here](https://docs.google.com/spreadsheets/d/1eu2Y2S0lyAUqxlvtnPBCO55OrxidYV7SwVRkqelPcKk/edit?gid=0#gid=0)

## Installation

Please use the following command to install this package in `R`

```
library(devtools)
devtools::install_github("satyanarayan-rao/SMTrackR")
```

## Generating Protein-DNA binding Visualization


### Using data from mouse 8Cell stage

The data is processed from [Wang et al., Nat Commun., 2021](https://pubmed.ncbi.nlm.nih.gov/33623021/) 

```
library(SMTrackR)
SMTrackR::plotFootprints(organism = "mmusculus", model = "8cell", condition = "WT", 
               genome_assembly = "mm10", type = "SMF", chromosome = "chr5", 
               start = "113847750",  stop  = "113847780", tr = "8cell", label = "remove_dup_true", 
               fp_cap = 50, remove_dup = T)
```
Please use the command below to generate a heatmap centered at `chr2L:480290-480320` in __Drosophila Melanogaster__ S2 cells.Here dSMF data is sourced from Krebs et al., 2017 Mol Cell.:

```
library(SMTrackR)
SMTrackR::plotFootprints()
```

