#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy
import gzip
from collections import defaultdict


# In[18]:


def loadFasta(filename):
    """ Parses a classically formatted and possibly 
        compressed FASTA file into a list of headers 
        and fragment sequences for each sequence contained.
        The resulting sequences are 0-indexed! """
    if (filename.endswith(".gz")):
        fp = gzip.open(filename, 'rb')
    else:
        fp = open(filename, 'rb')
    # split at headers
    data = fp.read().split('>')
    fp.close()
    # ignore whatever appears before the 1st header
    data.pop(0)     
    headers = []
    sequences = []
    for sequence in data:
        lines = sequence.split('\n')
        headers.append(lines.pop(0))
        sequences.append(''.join(lines))
    return (headers, sequences)


# In[19]:


filename1 = "GCA_009858895.3_ASM985889v3_genomic.fna.gz"
header1, seq1 = loadFasta(filename1)
genome = "+" + seq1[0]


# In[20]:


filename2 = "GCF_009858895.2_ASM985889v3_genomic.fna.gz"
header2,seq2 = loadFasta(filename2)


# In[21]:


header1


# In[22]:


header2


# In[23]:


seq1 == seq2


# In[24]:


def construct_nonoverlappingkmers(filename, k):
    def contig_name(headers):
        X = []
        for item in headers:
            X.append(item.split(' ')[0])
        return X

    def Construct_kmerprofile(X, sequences, k):
        '''Construct kmer profile given contig name and sequence
        parameter
        X - list of contig name
        sequences - list of reference sequences
        k - length of kmer
        '''
        # k = 45 # set the kmer length
        D = {}
        for i in range(len(sequences)): 
            contig = X[i]
            seq = sequences[i] # 0 offset
            num = len(seq)//k
            for i in range(num+1):
                D[contig] = D.get(contig, []) + [seq[i*k:(i+1)*k]]
        return D

    def main(filename, k):
        headers, sequences = loadFasta(filename)
        X = contig_name(headers)
        D = Construct_kmerprofile(X, sequences, k)
        return D
    
    D = main(filename, k)
    return D


# In[25]:


for k in range(5,12):
    kmers = construct_nonoverlappingkmers(filename1, k)
    print "kmerlength", k, "kmer num", len(kmers['MN908947.3'])


# In[26]:


for k in range(5,12):
    print k, len(seq1[0]) - k + 1, 4**k
print "minimal k is 8, maximal to be determined"


# In[27]:


# search kmer count
def mapping_position(genome, Kmer_Dict, k):   
    def create_kmer_profile(Kmer_Dict):
        KmerProfile = defaultdict(list)
        for sample,kmers in Kmer_Dict.iteritems():
            for seq in kmers:
                KmerProfile[seq]
                rev = ''.join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in reversed(seq)]) 
                KmerProfile[rev]
        print len(KmerProfile)
        return KmerProfile

    def mapping(seq, KmerProfile, k):
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i:i+k]

            if 'N' in kmer:
                continue

            c = KmerProfile.get(kmer, -1)

            if c < 0:
                continue

            KmerProfile[kmer] += [i]
        return KmerProfile

    KmerProfile = create_kmer_profile(Kmer_Dict)
    PositionProfile = mapping(genome, KmerProfile, k)
    return PositionProfile


# In[28]:


l = 10
PositionProfile = mapping_position(genome, kmers, k)


# In[29]:


# select anchor
# unique: forward = 1 reverse = 0. for all kmer list, check forward + reverse == 1
# not adjacent
def create_anchors(kmers, PositionProfile, k):
    def unique_kmers(kmers, PositionProfile):
        candidatelist = kmers.values()[0]
        candidates = []
        for kmer in candidatelist:
            assert PositionProfile[kmer] >0
            krev = ''.join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in reversed(kmer)]) 
            poslist = PositionProfile[kmer] + PositionProfile[krev]
            if len(poslist) == 1:
                candidates.append((kmer, PositionProfile[kmer][0]))
        return candidates

    #n =numpy.array(anchor_Dict['chr1'])
    def find_blocks(Indexlist, k):
        count = 0
        Blocks = []
        for i,s in enumerate(Indexlist):
            if i == 0:
                if Indexlist[i+1] - Indexlist[i] <= k:
                    count += 1
                    start = i
            elif i > 0 and i < len(Indexlist)-1:
                if Indexlist[i] - Indexlist[i-1] > k and Indexlist[i+1] - Indexlist[i] <= k:
                    count +=1
                    start = i
                elif Indexlist[i] - Indexlist[i-1] <= k and Indexlist[i+1] - Indexlist[i] > k:
                    end = i+1
                    Blocks.append((start, end))
            else:
                if Indexlist[i] - Indexlist[i-1] <= k:
                    end = i+1
                    Blocks.append((start, end))

        return count, Blocks


    def sort_by_kmerlength(poslist, k):
        poslist = numpy.array(poslist).astype(int)
        errorindex = []
        count, blocks = find_blocks(poslist, k)
        for s,e in blocks:
            if e - s < 3:
                for i in range(s+1,e):
                    errorindex.append(i)
            else:
                for i in range(s+1,e-1):
                    errorindex.append(i)
        return errorindex
    
    candidates = unique_kmers(kmers, PositionProfile)
    anchor_info = numpy.array(candidates)
    poslist = anchor_info[:,1]
    errorindex = sort_by_kmerlength(poslist, k)
    remaining = numpy.delete(anchor_info, errorindex, 0)
    return anchor_info, remaining


# In[30]:


anchors_info, remaining = create_anchors(kmers, PositionProfile, k)
numpy.save("anchorcandidates", anchors_info)


# In[67]:


len(anchors_info), len(remaining)


# In[68]:


poslist = anchor_info[:,1]


# In[69]:


errorindex = sort_by_kmerlength(poslist, k)


# In[71]:


errorindex
l = numpy.load("anchorcandidates.npy")


# In[72]:


len(l)


# In[ ]:




