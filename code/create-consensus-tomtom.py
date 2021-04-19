#!/usr/bin/env python
# coding: utf-8

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import argparse
import subprocess

TOMTOM='~/meme/bin/tomtom'
CHEN2MEME='/archive/gl/shared/software/meme-5.0.5/scripts/chen2meme'
MOUSEMEME='/archive/gl/shared/user/jhammelm/shared_data/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme'

# In[2]:
parser = argparse.ArgumentParser()
parser.add_argument('similarity')
parser.add_argument('motifs')
parser.add_argument('outfile')
parser.add_argument('-match','--match',action='store_true',default=False)
parser.add_argument('-cluster','--clustering',choices=['agglomerative','affinity'],default='affinity')
opts = parser.parse_args()

header = True
motifs = set()
for line in open(opts.similarity):
    if header:
        header = False
        continue
    if line[0] == '#':
        continue
    if line == '':
        break
    data = line.strip().split('\t')
    if data[0] != '':
        motifs.update([data[0]])


motifs = list(motifs)
similarity = np.zeros((len(motifs),len(motifs)))
offset = np.zeros((len(motifs),len(motifs)))
strand = -np.ones((len(motifs),len(motifs)))
header=True
for line in open(opts.similarity):
    if header:
        header = False
        continue
    if line[0] == '#':
        continue
    if line == '':
        break
    data = line.strip().split('\t')
    if data[0] != '':
        m1 = motifs.index(data[0])
        m2 = motifs.index(data[1])
        similarity[m1,m2] = -np.log10(float(data[5]))
        offset[m1,m2] = int(data[2])
        if data[9] == '+':
            strand[m1,m2]=1.0
            
similarity = np.nan_to_num(similarity)
similarity = np.clip(similarity,0,10)
similarity_cor = np.zeros(similarity.shape)
from scipy.stats import pearsonr
for i in range(similarity.shape[0]):
    for j in range(similarity.shape[1]):
        similarity_cor[i,j] = pearsonr(similarity[i,:],similarity[j,:])[0]


sns.clustermap(similarity_cor)
plt.savefig(opts.outfile+'-correlation.svg')


# In[67]:


def align(cluster,offsets,similarity,strand):
    alignments = np.zeros((len(cluster),),dtype='int')
    fw = np.zeros((len(cluster),))
    first = sorted(cluster,key = lambda x: np.max(similarity[x,cluster[cluster != x]]))[0]
    total = set([first])
    alignments[cluster==first] = 0
    fw[cluster==first]=1
    while len(total) != len(cluster):
        remaining = [c for c in cluster if c not in total]
        current = sorted(remaining,key = lambda x: np.max(similarity[x,np.array(list(total),dtype='int')]))[-1]
        other = np.array([c for c in total if c != current],dtype='int')
        closest = other[np.argmax(similarity[current,other])]
        alignments[cluster==current] = alignments[cluster==closest] - offset[current,closest]
        fw[cluster==current] = strand[current,closest]*fw[cluster==closest]
        total.update([current])

    return alignments,fw


# In[59]:


def reverse_complement(string):
    rc = ""
    alpha = {'A':'T','C':'G','T':'A','G':'C'}
    for c in string:
        rc += alpha[c]
    return rc[::-1]


def motif_entropy(pcm):
    entropy=0
    for p in range(pcm.shape[0]):
        pwm = pcm[p,:]/np.sum(pcm[p,:])
        for c in range(pcm.shape[1]):
            if pwm[c] != 0:
                entropy += -pwm[c]*np.log2(pwm[c])
    return entropy/pcm.shape[0]

motif_dict = {}
for motif in open(opts.motifs).read().split('>')[1:]:
    #print(motif)
    motif_lines = motif.strip().split('\n')
    header = motif_lines[0].split()[0]
    seq_mat = np.zeros((len(motif_lines)-1,4))
    for i in range(1,len(motif_lines)):
        seq_mat[i-1,:] = np.array([float(v) for v in motif_lines[i].split('\t')],dtype='float')
    motif_dict[header] = seq_mat

distance = 1.0 - np.clip(similarity_cor,0,np.max(similarity_cor))
from sklearn.cluster import AgglomerativeClustering,AffinityPropagation
if opts.clustering == 'agglomerative':
    for thresh in [0.2,0.5,0.8]:
        agg_clust = AgglomerativeClustering(n_clusters=None,
                                            distance_threshold=1.0-thresh,
                                            affinity='precomputed',
                                            linkage='complete').fit(distance)
        seq_motifs = {}
        print('# clusters:',len(set(agg_clust.labels_)))
        for l in set(agg_clust.labels_):
            if np.where(agg_clust.labels_==l)[0].shape[0] == 0:
                continue
            if np.where(agg_clust.labels_==l)[0].shape[0] == 1:
                mind = np.where(agg_clust.labels_==l)[0][0]
                seq_motifs['m'+str(l+1)] = motif_dict[motifs[mind]]
            else:
                algns,fw = align(np.where(agg_clust.labels_==l)[0],offset,similarity,strand)
                max_left_shift = np.min(algns)
                algns_left = algns - max_left_shift
                motif_lens = [motif_dict[motifs[motif_i]].shape[0] for motif_i in np.where(agg_clust.labels_==l)[0]]
                all_seqs = np.zeros((np.max(algns_left)+max(motif_lens),4))
            
                for within_i,motif_i in enumerate(np.where(agg_clust.labels_==l)[0]):
                    if fw[within_i] == -1:
                        ro = [3,2,1,0]
                        motif_version = motif_dict[motifs[motif_i]][::-1,ro]
                    else:
                        motif_version = motif_dict[motifs[motif_i]]
                    for mpos,pos in enumerate(range(algns_left[within_i],
                                                    algns_left[within_i]+motif_dict[motifs[motif_i]].shape[0])):
                        all_seqs[pos,:] += motif_version[mpos,:]
                seq_motifs['m'+str(l+1)] = all_seqs
            
        with open(opts.outfile+'-aggthresh='+str(thresh)+'.motifs','w') as f:
            for key in seq_motifs.keys():
                f.write('>'+key +'\n')
                for i in range(seq_motifs[key].shape[0]):
                    f.write('\t'.join([str(v) for v in seq_motifs[key][i,:]])+'\n')
        subprocess.call([CHEN2MEME+' < '+opts.outfile+'-aggthresh='+str(thresh)+'.motifs > '+opts.outfile+'-aggthresh='+str(thresh)+'.meme'],shell=True)
        if opts.match:
            subprocess.call([TOMTOM+' -oc '+opts.outfile + '-aggthresh='+str(thresh)+' '+opts.outfile+'-aggthresh='+str(thresh)+'.meme '+MOUSEMEME],shell=True)
else:
    #opts.clustering == 'representative'
    #shrink motif entropy to end up with fewer clusters?
    preference = np.array([1.0/motif_entropy(motif_dict[motif]) for motif in motifs])*0.25
    
    agg_clust = AffinityPropagation(affinity='precomputed',preference=preference).fit(np.clip(similarity_cor,0,np.max(similarity_cor)))
    seq_motifs = {}
    print('# clusters:',len(set(agg_clust.labels_)))
    for label_ind,cluster_center in enumerate(agg_clust.cluster_centers_indices_):
        seq_motifs[motifs[cluster_center]] = motif_dict[motifs[cluster_center]]
        with open(opts.outfile+'-'+motifs[cluster_center]+'.cluster','w') as f:
            for mind in np.where(agg_clust.labels_==label_ind)[0]:
                f.write(motifs[mind]+'\n')
    with open(opts.outfile+'-affinityprop.motifs','w') as f:
        for key in seq_motifs.keys():
            f.write('>'+key +'\n')
            for i in range(seq_motifs[key].shape[0]):
                f.write('\t'.join([str(v) for v in seq_motifs[key][i,:]])+'\n')
    subprocess.call([CHEN2MEME+' < '+opts.outfile+'-affinityprop.motifs > '+opts.outfile+'-affinityprop.meme'],shell=True)
    if opts.match:
        subprocess.call([TOMTOM+' -oc '+opts.outfile + '-affinityprop '+opts.outfile+'-affinityprop.meme '+MOUSEMEME],shell=True)
    
            
