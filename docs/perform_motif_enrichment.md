## Mouse Motif Enrichment
We created a [script](https://raw.githubusercontent.com/gifford-lab/ReprogrammingRecovery/main/perform_motif_enrichment/perform_motif_enrichment.py) that takes as input on MACS2 narrowPeak formatted bed files with mm10 mouse genome coordinates and performs processing of input regions to select regions, selection negative (control sequences, and runs motif enrichment using AME, DREME, KMAC, and HOMER. 

This script has been tested for Ubuntu v18.04.3 with python v3.7.9. To use the script, you will have to modify the variables to the correct locations of mm10_fasta, ame_path, meme_path, kmac_path, and homer_path within your compute environment.

![](/assets/motif_enrichment.png)

### Installation
First install data and software dependencies:
- [mouse genome mm10](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz)
- [MEME Suite 5.0.5](https://meme-suite.org/meme/doc/download.html) 
- [GEM 3.4](http://groups.csail.mit.edu/cgs/gem/kmac/)
- [HOMER 4.9.1](http://homer.ucsd.edu/homer/)
- [bedtools 2.29.2](https://bedtools.readthedocs.io/en/latest/content/installation.html)

Then to perform motif enrichment, download the python script and related files:
```
git clone https://github.com/gifford-lab/ReprogrammingRecovery.git
cd ReprogrammingRecovery
```
In the perform_motif_enrichment folder, open the python script perform_motif_enrichment.py and edit paths to the paths to data and software dependencies within your local compute environment.

To run the motif enrichment script, you must also bedtools to your path (i.e. export PATH=/path/to/my/bin/bedtools:$PATH).

Install should be complete. To see program arguments: 
```
cd perform_motif_enrichment
python perform_motif_enrichment.py -h
```

## Run on example data
We will show an example of how to use our script to run all four motif discovery methods on accessibility data from dopaminergic midbrain neurons. Download example [MACS2 peak file](http://reprogramdata.csail.mit.edu/atac/peaks/DP_MB-dopaminergic_rep1-threshold0.01.filt.narrowPeak.gz). Then to perform motif enrichment on the top 500 peaks with multitissue enhancer background sequences:
```
python perform_motif_enrichment.py DP_MB-dopaminergic_rep1-threshold0.01.filt.narrowPeak.gz -db all -bg ../data/mm10_shared_enhancers.fa -n 500 -ame -homer -meme -kmac -o MB-dopaminergic
```
Output will be four folders representing results from AME, MEME (DREME), KMAC, and HOMER. Each methods output folder contains an html file which summarizes the results. 

Runtime should be less than 1 hour on a normal desktop computer. 


