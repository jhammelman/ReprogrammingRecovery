# Reprogramming Recovery

![](docs/assets/reprogramming_schematic.png)

## This git repo contains code and data and the [website](https://cgs.csail.mit.edu/ReprogrammingRecovery/) supporting the paper "Ranking Reprogramming Factors for Directed Differentiation."

## Citation
[Ranking Reprogramming Factors for Directed Differentiation](https://www.biorxiv.org/content/10.1101/2021.05.14.444080v2)\
Jennifer Hammelman, Tulsi Patel, Michael Closser, Hynek Wichterle, David Gifford\
bioRxiv 2021.05.14.444080; doi: https://doi.org/10.1101/2021.05.14.444080 

## cluster_consensus_motifs
Folder contains methods for generating a consensus transcription factor motif data base of 107 consensus motifs representing 356 motifs in the HOCOMOCOv11 core mouse database using the script 00_cluster_motifs.sh.

## perform_motif_enrichment
Folder contains scripts for performing motif enrichment using DREME, AME, HOMER, and KMAC. 

## data
### consensus_HOCOMOCOv11_core_MOUSE-affinityprop-handannotated.meme
consensus database of 107 motifs representing TF families in meme format

### consensus_HOCOMOCOv11_core_MOUSE-affinityprop-handannotated.motifs
consensus database of 107 motifs representing TF families in chan PWM format

### consensus_HOCOMOCOv11_core_MOUSE-affinityprop-mapping.txt
mapping of 356 transcription factors to 107 motif consensus families

### mouse_ensemble_tfs_from_lambertetal_isyes.unique.txt
list of ensemble mouse gene names of 1,374 transcription factors curated in Lambert et al., 2018.
