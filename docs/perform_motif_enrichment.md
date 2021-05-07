## Mouse Motif Enrichment
We created a [script](https://raw.githubusercontent.com/gifford-lab/ReprogrammingRecovery/main/perform_motif_enrichment/perform_motif_enrichment.py) that takes as input on MACS2 narrowPeak formatted bed files with mm10 mouse genome coordinates and performs processing of input regions to select regions, selection negative (control sequences, and runs motif enrichment using AME, DREME, KMAC, and HOMER. 

This script has been tested for Ubuntu v18.04.3 with python v3.7.9. To use the script, you will have to modify the variables to the correct locations of mm10_fasta, ame_path, meme_path, kmac_path, and homer_path within your compute environment.

![](/assets/motif_enrichment.png)


