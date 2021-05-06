## Consensus Mouse Motif Database
We generated and hand-annotated a consensus motif database of 107 mouse transcription factor groups by clustering the HOCOMOCOv11 mouse database of 356 transcription factor motifs.

![](/assets/consensus_database_method.png)

The script to compute similarity and perform affinity propagation clustering of the motifs is available [here](https://github.com/jhammelman/ReprogrammingRecovery/blob/main/cluster_consensus_motifs/00_cluster_motifs.sh).

We provide this database in the following formats: 
1. [MEME](https://raw.githubusercontent.com/jhammelman/ReprogrammingRecovery/main/data/consensus_HOCOMOCOv11_core_MOUSE-affinityprop-handannotated.meme)
2. [Chan format PCM](https://raw.githubusercontent.com/jhammelman/ReprogrammingRecovery/main/data/consensus_HOCOMOCOv11_core_MOUSE-affinityprop-handannotated.motifs)

We also provide a [mapping](https://raw.githubusercontent.com/jhammelman/ReprogrammingRecovery/main/data/consensus_HOCOMOCOv11_core_MOUSE-affinityprop-renamed-cluster-mapping.txt) from consensus motif to transcription factor motifs within the cluster. Stars denote the transcription factor motif within the group that is used as the representative motif for the cluster. 
