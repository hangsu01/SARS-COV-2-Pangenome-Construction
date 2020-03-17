#!/usr/bin/env python
# coding: utf-8

# In[6]:


import CCGG_extension as CCGG
import sys
import numpy
import gzip
from collections import defaultdict
import pandas as pd
import json
import gzip
from tqdm import tqdm_notebook
#sys.path.append('/nas/longleaf/home/hangsu/.local/lib/python3.7/site-packages')

# import sys
# !{sys.executable} -m pip install python-levenshtein

import Levenshtein
#from ipywidgets import IntProgress
#import PIL
#import Levenshtein


# In[7]:


#sys.path.append('/nas/longleaf/home/hangsu/.local/lib/python3.7/site-packages')
# import CCGG_extension as CCGG
# sys.prefix


# In[8]:


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

def writeGraphFasta(filename, input_dict, keylist=["src", "dst"]):
    """Write the given node or edge file as a FASTA file. Overwrites existing file. Will create file if it doesn't exist. 
    def writeFasta(self:<GraphicalGenome>, filename:<str>, input_dict:<dict>, keylist=["src","dst"]:<list[str]>) -> Void

    Parameters:
        filename: <str> - absolute path for the file you wish to write. 
        input_dict: <dict> - Node or Edge dictionary you wish to write. 
        keylist: <list[str]> - list of strings of dictionary keys that you want to ignore during write. 
                This will require you to write these values on your own or just discard them. 
    """
    sorted_keys = sorted(input_dict.keys()) 
    with open(filename, "w+") as fastafile:
        # If iterating through the edges, write the edges in the correctly ordered format
        if (sorted_keys[0][0] == "E"):
            for edge in sorted_keys:
                # If header has not been evaluated, just re-write the header wholesale without any analysis
                if "hdr" in input_dict[edge].keys():
                    line = ">" + edge + ";" + input_dict[edge]["hdr"] + "\n"
                    line += input_dict[edge]["seq"] + "\n"
                    continue
                line = ">" + edge + ";{" 
                # Source
                line += '"src":"' + input_dict[edge]["src"] + '",'
                # Destination
                line += '"dst":"' + input_dict[edge]["dst"] + '"'
                for key in input_dict[edge].keys():
                    if key == "seq":
                        continue
                    if key in keylist:
                        continue
                    line += ',"' + key + '":' + json.dumps(input_dict[edge][key], separators=(",", ":"))
                line += "}\n"
                line += input_dict[edge]["seq"] + "\n"
                fastafile.write(line)
        # If iterating over nodes, just write the nodes normally
        else:
            for i in sorted_keys:
                line = ">" + i + ";"
                obj = {}
                for j in input_dict[i].keys():
                    if j == 'seq':
                        continue
                    obj[j] = input_dict[i][j]
                line += json.dumps(obj, separators=(",", ":"))
                line += "\n" + input_dict[i]['seq'] + "\n"
                fastafile.write(line)


# In[9]:


# search kmer count
def mapping_position(genome, anchor_filename, k):   
    def create_kmer_profile(kmerlist):
        KmerProfile = defaultdict(list)
    #         for sample,kmers in Kmer_Dict.iteritems():
        for seq in kmerlist:
            KmerProfile[seq]
            rev = ''.join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in reversed(seq)]) 
            KmerProfile[rev]
    #         print len(KmerProfile)
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

    # determine uniqueness
    def unique_kmers(kmerlist, PositionProfile):
        candidates = []
        for kmer in kmerlist:
            if len(PositionProfile[kmer]) > 0:
                krev = ''.join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in reversed(kmer)]) 
                poslist = PositionProfile[kmer] + PositionProfile[krev]

                if len(poslist) == 1:
                    #print PositionProfile[kmer]
                    candidates.append((kmer, PositionProfile[kmer][0]))
        return candidates

    # determine monotonicity
    def binary_search(arr, val, l, r): 
        if l == r: 
            if arr[l] > val: 
                return l 
            else: 
                return l+1
        if l > r: 
            return l 

        mid = (l+r)/2
        if arr[mid] < val: 
            return binary_search(arr, val, mid+1, r) 
        elif arr[mid] > val: 
            return binary_search(arr, val, l, mid-1) 
        else: 
            return mid 

    # NlogN 
    def efficientDeletionSort(array):
        subsequence_end = [0] # index of the element
        predecessors = [-1] # predecessor index
        for i in range(1,len(array)):
            arr = array[i]
            # can do binary search instead, just skip to make it faster
            if arr > array[subsequence_end[-1]]:
                predecessors += [subsequence_end[-1]]
                subsequence_end += [i]
            else:
                # preform binary search
                minimum_end = [array[j] for j in subsequence_end] # element in current subsequence

                insert_point = binary_search(minimum_end, arr, 0, len(minimum_end)-1)
                if insert_point > 0:
                    predecessors += [subsequence_end[insert_point-1]]
                else:
                    predecessors += [-1]

                if insert_point > len(subsequence_end)-1: # arr == array[subsequence_end[-1]] 
                    subsequence_end += [i]
                elif arr < array[subsequence_end[insert_point]]:
                    subsequence_end[insert_point] = i 

        # backtrack
        pre = subsequence_end[-1]
        listIndex = []
        while pre != -1:
            listIndex.append(pre)
            pre = predecessors[pre]
        listIndex.reverse()
        longest_subsequence = [array[i] for i in listIndex]
        return listIndex, longest_subsequence

    def monotonic_kmers(candidates):
        candidates = numpy.array(candidates)
        poslist = candidates[:,1].astype(int)
        listIndex, longest_subsequence = efficientDeletionSort(poslist)
        monotonic_anchor = candidates[listIndex,:]
        return monotonic_anchor

    anchor_c = numpy.load(anchor_filename)
    kmerlist = anchor_c[:,0]
    KmerProfile = create_kmer_profile(kmerlist)
    PositionProfile = mapping(genome, KmerProfile, k)
    candidates = unique_kmers(kmerlist, PositionProfile)
    final = monotonic_kmers(candidates)
    
    return final


# In[10]:


def getsequence_info(genomefile):
    def contig_name(headers):
        X = []
        for item in headers:
            X.append(item.split('|')[0])
        return X
    def contig_name_CNCB(headers):
        X = []
        for item in headers:
            X.append(item.split('|')[-1][1:])
        return X
    
    header, sequence = loadFasta(genomefile)
    if genomefile.startswith('CNCB'):
        samples = contig_name_CNCB(header)
    else:
        samples = contig_name(header)
    
    for i in range(len(samples)):
        if samples[i] == 'GWHABKP00000000':
            print header[i]
            del samples[i]
            del header[i]
            del sequence[i]
            break   
    return header, sequence, samples

def integrate_info(genomefile, anchorfile, k):
    header, sequence, samples = getsequence_info(genomefile)
    anchor_c = numpy.load(anchorfile)
    Anchor_Info = {}
    for i in range(len(header)):
        genome = "+" + sequence[i].upper()
        # if not complete genome
        if len(genome) < 10000:
            continue
        # if too many ambiguous bases
        basenum = genome.count("A") + genome.count("G") + genome.count('C') + genome.count("T")
        if float(basenum)/len(genome) < 0.9:
            continue
        anchormapping = mapping_position(genome, anchorfile, k)

        Final_Dict= dict(anchormapping) # current mapping

        samplename = samples[i]
        D = {}
        for anchorseq, refpos in anchor_c:
            anchorname = "A%05d" % (int(refpos)/k + 1)
            D[anchorname] = Final_Dict.get(anchorseq, "?")
        Anchor_Info[samplename] = D
        
    return Anchor_Info
    


# In[46]:


genomefile = "COVID_19sequences.fasta"
genomefile = 'Covid19sequences0304.fasta'

anchorfile = "anchorcandidates_11.npy"
k = 11
anchor_info = integrate_info(genomefile, anchorfile, k)
header, sequence = loadFasta(genomefile)


# In[47]:


#print len(anchor_info)

def dynamic_construct(anchor_info, strainlist = []):
    df = pd.DataFrame(anchor_info)
    if len(strainlist) == 0:
        d = df == '?'
        Anchorlist = df[d.sum(axis = 1) == 0]
        return Anchorlist
    else:
        df = df.loc[:,strainlist]
        d = df == '?'
        Anchorlist = df[d.sum(axis = 1) == 0]
        return Anchorlist


# In[48]:


Anchorlist = dynamic_construct(anchor_info)


# In[ ]:





# In[50]:


# collapse adjacent anchors
def collapsed_adjacent(AnchorDataFrame, k):
    
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

    columnnames = AnchorDataFrame.columns
    for i in columnnames:
        anchors = AnchorDataFrame.index
        poslist = AnchorDataFrame[i].values.astype(int)
        assert sum(sorted(poslist) == poslist) == len(poslist)
        errorindex = sort_by_kmerlength(poslist, k)
        index = [i for i in range(len(anchors)) if i not in errorindex]
        AnchorDataFrame = AnchorDataFrame.iloc[index,:]
        
    return AnchorDataFrame


# In[51]:


def create_pangenome(AnchorDataFrame, genomefile,k):
    AnchorDataFrame = collapsed_adjacent(AnchorDataFrame,k)
    
    def genome_info(genomefile):
        header, sequence, samples = getsequence_info(genomefile)
        genome_Dict = {}
        for i in range(len(header)):
            genome = "+" + sequence[i]
            genome_Dict[samples[i]] = genome
        return genome_Dict
    
    def create_nodefile(genome_Dict, AnchorDataFrame,k):
        genome_Dict = genome_info(genomefile)
        Anchors = {}
        refname = 'MN908947'
        for anchor in AnchorDataFrame.index:
            genome = genome_Dict[refname]
            Anchors[anchor] = dict(AnchorDataFrame.loc[anchor])
            pos = Anchors[anchor][refname]
            Anchors[anchor]['seq'] = genome[int(pos):int(pos) + k]  
        return Anchors
    
    def get_all_edge_sequence(AnchorDataFrame, genome_Dict, src, dst):
        Strainlist = AnchorDataFrame.columns
        D = {}
        for strain in Strainlist:
            seq = genome_Dict[strain].upper()
            if src == 'SOURCE':
                for i, s in enumerate(seq):
                    if s!= "N" and s != '+':
                        start = i
                        break
            else:
                start = int(AnchorDataFrame.loc[src, strain]) + k

            if dst == 'SINK':
                end = len(seq)
            else:
                end = int(AnchorDataFrame.loc[dst, strain])

            edge_seq = seq[start:end]
            D[strain] = edge_seq

        return D
    
    def get_edge_info(src,dst, AnchorDataFrame, genome_Dict):
        D = get_all_edge_sequence(AnchorDataFrame, genome_Dict, src, dst)
        Strain = Anchorlist.columns
        Edge_D = {}
        index = 0
        my_k = 0
        while my_k in range(len(Strain)):
            index += 1
            if src != "SOURCE":
                edgename = 'E%06d' % (int(src[1:])*100 + index)
            else:
                edgename = 'S%06d' % (index)
            Edge_D[edgename] = {}
            strain = Strain[my_k] 
            Name = []        
            Name.append(strain)

            if len(Strain) - my_k>1:
                for j in range(my_k+1, len(Strain)):
                    strain1 = Strain[j]
                    if D[strain] == D[strain1]:
                        Name.append(strain1)
            Strain = [x for x in Strain if x not in Name]

            Edge_D[edgename]['seq'] = D[strain]
            Edge_D[edgename]['strain'] = sorted(Name)
            Edge_D[edgename]['src'] = src
            Edge_D[edgename]['dst'] = dst

        return Edge_D    
    
    def create_edgefile(genome_Dict, AnchorDataFrame, k):
        nodelist = ["SOURCE"] + list(AnchorDataFrame.index) + ['SINK']
        Edge_info = {}
        total = 0
        for i in range(len(nodelist) - 1):
            src = nodelist[i]
            dst = nodelist[i + 1]
            Edges = get_edge_info(src,dst, AnchorDataFrame, genome_Dict)
            total += len(Edges)
            Edge_info.update(Edges)    
        return Edge_info
    
    genome_Dict = genome_info(genomefile)
    Anchors = create_nodefile(genome_Dict,AnchorDataFrame,k)
    Edges = create_edgefile(genome_Dict, AnchorDataFrame, k)
    
    return Anchors, Edges


# In[52]:


Anchor_Dict,Edge_Dict = create_pangenome(Anchorlist, genomefile, 11)
writeGraphFasta("Cov_Nodefile.fa",Anchor_Dict)
writeGraphFasta("Cov_Edgefile.fa",Edge_Dict)


# In[53]:


import CCGG_extension as CCGG


# In[54]:


graph = CCGG.GraphicalGenome("Cov_Nodefile.fa", "Cov_Edgefile.fa", nnamelen= 6, enamelen=7)


# # add annotation

# In[69]:


# add genome feature
# [seqid, source, type, start, end, score, strand, phase, attribute]
def get_gff_data(filename):
    with gzip.open(filename, 'rb') as fp:
        data = fp.readlines()
    return data

def parse_feature(data):
    Parser = {}
    count = 0
    for line in data:
        values = line.decode().split('\t')
        if len(values) < 9:
            continue
        start, end = int(values[3]), int(values[4])
        # parse attributes
        attributes = values[8]
        Att = {}
        attributes_info = attributes.split(';')
        for item in attributes_info:
            k, v = item.split('=')
            Att[k] = v
        #assert "ID" in Att.keys()
        # build in   
        name = Att["ID"]
        Parser[name] = {}
        Parser[name]["start"] = start
        Parser[name]['end'] = end
        Parser[name]['att'] = Att
        Parser[name]['strand'] = values[6]
        Parser[name]['type'] = values[2]
    return Parser

def annotatefeature(graph, filename,k):    
    def write_cigar(g_start, g_end, i_start, i_end):
        ''' input genestart, geneend, itemstart, item_end Gene include both ends [], items include the start [) '''
        g_start = int(g_start)
        g_end = int(g_end)
        i_start = int(i_start)
        i_end = int(i_end)
        if i_start >= g_start and i_end - 1 <= g_end:
            Cigar = str(i_end - i_start)+'M'

        elif i_start < g_start and i_end - 1 <= g_end:
            miss = g_start - i_start
            Cigar = str(miss)+'S'+str(i_end-1-g_start+1)+'M'

        elif i_start < g_start and i_end -1 > g_end: # gene inside a Node
            miss = g_start - i_start
            Cigar = str(miss)+'S'+str(g_end - g_start + 1)+'M' + str(i_end - 1 - g_end)+'S'

        elif i_end - 1 > g_end:
            miss = i_end -1 - g_end
            Cigar = str(g_end - i_start + 1)+'M'+str(miss)+'S'

        return Cigar


    def find_itempos(Graph, item, k):
        '''Input graph item, return start and end pos [istart, iend), end pos not included'''
        if item.startswith('A'):
            i_start = int(Graph.nodes[item]['MN908947'])
            i_end = i_start + k
        elif item.startswith('F') or item.startswith('B') or item == 'SOURCE' or item == 'SINK':
            return 
        else:
            snode = Graph.edges[item]['src']
            if snode == 'SOURCE':
                i_start = 1
            else:
                i_start = int(Graph.nodes[snode]['MN908947']) + k

            i_end = i_start + len(Graph.edges[item]['seq'])

        return (i_start, i_end)

    def addannotation(Graph, Feature_info,k):
        featurelist = Feature_info.keys()
        for g_ID in featurelist:
            D = Feature_info[g_ID]
            g_start = int(D['start'])
            g_end = int(D['end'])
            strand = D['strand']
            key = D['type']
            sanchor, eanchor = Graph.boundingAnchors('MN908947', g_start, g_end)
            itemlist = Graph.tracePath('MN908947', sanchor, eanchor)

            for item in itemlist:
                if item == 'SOURCE' or item == 'SINK':
                    continue
                i_start, i_end = find_itempos(graph, item,k)
                if i_end <= g_start or i_start > g_end:
                    continue
                    
                cigar = write_cigar(g_start, g_end, i_start, i_end)
#                 if cigar == '0M':
#                     print( cigar, item,g_start, g_end, i_start, i_end, len(Graph.edges[item]['seq']))
                
                if item.startswith('A'):
                    annote = '|'.join((g_ID, strand, cigar))
                    if annote in Graph.nodes[item].get(key, []):
                        continue
                    Graph.nodes[item][key] = Graph.nodes[item].get(key, []) + ['|'.join((g_ID, strand, cigar))]
                else:
                    annote = '|'.join((g_ID, strand, cigar))
                    if annote in Graph.edges[item].get(key, []):
                        continue
                    Graph.edges[item][key] = Graph.edges[item].get(key, []) + ['|'.join((g_ID, strand, cigar))]
        return Graph
    
    data = get_gff_data(filename)
    feature_info = parse_feature(data)
    graph = addannotation(graph, feature_info,k)
    return graph

annotationfile1 = "GCF_009858895.2_ASM985889v3_genomic.gff.gz"
# data1 = get_gff_data(annotationfile1)
# annotationfile2 = "GCA_009858895.3_ASM985889v3_genomic.gff.gz"
# data2 = get_gff_data(annotationfile2)
graph = annotatefeature(graph, annotationfile1,11)


# In[131]:


writeGraphFasta("Cov_Nodefile.fa",graph.nodes)
writeGraphFasta("Cov_Edgefile.fa",graph.edges)


# # add variants

# In[11]:


#pip install --user python-Levenshtein


# In[71]:


import Levenshtein
def addvariants(Graph, maxlength,  refstrain = "MN908947"):
    def makeCigar(seq, ref):
        if (len(seq) > 16384) or (len(ref) > 16384):
            rmid = len(ref)/2        
            smid = len(seq)/2
            prox = makeCigar(seq[:smid],ref[:rmid])
            dist = makeCigar(seq[smid:],ref[rmid:])
            return prox+dist
        ops = Levenshtein.editops(seq,ref)
        code = ['=' for i in xrange(len(seq))]
        offset = 0
        for op, si, di in ops:
            if (op == "replace"):
                code[si+offset] = 'X'
            elif (op == "insert"):
                code.insert(si+offset,'D')
                offset += 1
            elif (op == "delete"):
                code[si+offset] = 'I'# LM: fixed bug here 2019-04-15
        cigar = ''
        count = 1
        prev = code[0]
        for c in code[1:]:
            if (c == prev):
                count += 1
            else:
                cigar += "%d%c" % (count, prev)
                count = 1
                prev = c
        cigar += "%d%c" % (count, prev)
        return cigar

    # check the length of the cigar string is valid
    def search_valid_size(variants,seq):
        ref_length,alt_length = 0,0
        for v in variants:
            base, op = v
            base = int(base)
            if op == '=' or op == 'X':
                ref_length += base
                alt_length += base
            elif op == 'I':
                alt_length += base
            elif op == 'D':
                ref_length += base

        return len(seq) == alt_length 

    # check the length of the cigar string is valid
    def valid_size(variants,seq,reference):
        ref_length,alt_length = 0,0
        for v in variants:
            base, op = v
            base = int(base)
            if op == '=' or op == 'X':
                ref_length += base
                alt_length += base
            elif op == 'I':
                alt_length += base
            elif op == 'D':
                ref_length += base

        assert len(seq) == alt_length 
        assert len(reference) == ref_length

    def valid_eqaul_and_mismatch(variants,seq,reference):
        ref_pos,alt_pos = 0,0
        for v in variants:
            base,op = v
            base = int(base)
            if op == '=':
                assert seq[alt_pos:alt_pos+base] == reference[ref_pos:ref_pos+base]
                ref_pos += base
                alt_pos += base
            elif op == 'X':
                assert seq[alt_pos:alt_pos+base] != reference[ref_pos:ref_pos+base]
                ref_pos += base
                alt_pos += base
            elif op == 'I':
                alt_pos += base
            elif op == 'D':
                ref_pos += base

    #assert len(reference) == ref_length
    def split(variant):
        splitVariants = []
        previous = 0
        for i in range(len(variant)):
            v = variant[i]
            if v.isdigit():
                continue
            else:
                splitVariants += [variant[previous:i]]
                splitVariants += [v]
                previous = i + 1

        numbers, types = [],[]
        for i in range(len(splitVariants)):
            v = splitVariants[i]
            if i %2 == 0:
                numbers += [v]
            else:
                types += [v]

        variant = []
        for i in range(len(numbers)):
            variant += [(numbers[i],types[i])]
        return variant

    def findBedge(Graph, src, dst, refstrain):
        edgelist = Graph.outgoing[src]
        for edge in edgelist:
            sequence = Graph.edges[edge]['seq'].upper()
            strain = Graph.edges[edge]['strain']
            if refstrain in strain:
                B_dst = Graph.edges[edge]['dst']
                if B_dst == dst:
                    return edge
                else:
                    return ""
        else:
            return ''

    # mask Ns on the alternative path
    def maskNs(seq):
        seq = seq.upper()
        newseq = ''
        for i in seq:
            if i not in "AGCT":
                newseq += i.lower()
            else:
                newseq += i
        return newseq
    
    def add_variants(Graph, maxlength, refstrain, gapopen = False):
        Nodelist = sorted(Graph.nodes.keys())
        Nodelist = ["SOURCE"] + Nodelist
        for node in Nodelist:
            sanchor = node
            eanchor = Graph.nextAnchor(sanchor)
            edgelist = Graph.outgoing[sanchor]
            ref_edge = findBedge(Graph, sanchor, eanchor, refstrain)

            # fix redundant Ns on SOURCE
            if sanchor == "SOURCE":
                seq = Graph.edges[ref_edge]['seq'].upper()
                for i, s in enumerate(seq):
                    if s!= "N":
                        start = i
                        break
                ref_seq = seq[start:]
                ref_length = len(ref_seq)
                Graph.edges[ref_edge]['seq'] = ref_seq  
            # fix redundant Ns on SINK
            elif eanchor == "SINK":
                seq = [i for i in reversed(Graph.edges[ref_edge]['seq'].upper())]
                for i, s in enumerate(seq):
                    if s!= "N":
                        start = i
                        break         
                ref_seq = seq[start:]
                ref_seq = ''.join([i for i in reversed(ref_seq)])
                ref_length = len(ref_seq)
                Graph.edges[ref_edge]['seq'] =  ref_seq
            else:
                ref_seq = Graph.edges[ref_edge]['seq'].upper()
                ref_length = len(ref_seq)



            if ref_length > maxlength:
                continue

            for alt_edge in edgelist:
                if alt_edge == ref_edge:
                    Graph.edges[ref_edge]['variants'] = str(ref_length) + '='
                    continue

                alt_seq = maskNs(Graph.edges[alt_edge]['seq']) # mask Ns
                alt_length = len(alt_seq)

                if alt_length > maxlength:
                    continue

                if gapopen:
                    delta = abs(alt_length - ref_length)
                    if delta > 100:
                        if delta > alt_length or delta > ref_length:
                            cigar = GapOpenAligner(ref_seq, alt_seq)
                            variants = split(cigar)
                            valid_size(variants, alt_seq, ref_seq)
                            valid_eqaul_and_mismatch(variants, alt_seq, ref_seq)
                        else:
                            cigar = makeCigar(alt_seq, ref_seq)
                            variants = split(cigar)
                            valid_size(variants, alt_seq, ref_seq)
                            valid_eqaul_and_mismatch(variants, alt_seq, ref_seq)
                    else:
                        cigar = makeCigar(alt_seq, ref_seq)
                        variants = split(cigar)
                        valid_size(variants, alt_seq, ref_seq)
                        valid_eqaul_and_mismatch(variants, alt_seq, ref_seq)
                else:
                    cigar = makeCigar(alt_seq, ref_seq)
                    variants = split(cigar)
                    valid_size(variants, alt_seq, ref_seq)
                    valid_eqaul_and_mismatch(variants, alt_seq, ref_seq)

                Graph.edges[alt_edge]['variants'] = cigar
    
    add_variants(Graph, maxlength, refstrain, gapopen = False)
    return Graph


# In[134]:


graph = addvariants(graph, 3000)
writeGraphFasta("Cov_Nodefile.fa",graph.nodes)
writeGraphFasta("Cov_Edgefile.fa",graph.edges)


# In[61]:


founder = []
for edge in edgelist:
    founder += graph.edges[edge]['strain']
    print graph.edges[edge]


# In[18]:


print set(founder)
#graph.getSubPath("SOURCE", 'A00005', init_strain= set(founder))
print 'A00005' in graph.nodes.keys()


# # validation

# In[63]:


def Validate_sequence(Graph, genome_Dict, strain, path = 0):
    letter2strain = {'A':'AJ', 'B':'B6', 'C':'129S1', 'I':'DBA2', 'D':"NOD", 'E':"NZO", 'F':"CAST", 'G':"PWK", 'H':"WSB"}
    conseq = Graph.reconstructSequence(strain, path)
#     with open('/csbiodataxw/KeaneGenomes/genomes/%s/Chr%s.seq' % (letter2strain[strain],str(chromo))) as fp:
    lineargenome = genome_Dict[strain][1:].upper()
    return conseq.upper() == lineargenome.upper()

strainlist = Anchorlist.columns
print len(strainlist)
def genome_info(genomefile):
    header, sequence, samples = getsequence_info(genomefile)
    genome_Dict = {}
    for i in range(len(header)):
        genome = "+" + sequence[i]
        genome_Dict[samples[i]] = genome
    return genome_Dict

genome_Dict = genome_info(genomefile)
for strain in strainlist:
    print Validate_sequence(graph, genome_Dict, strain)


# In[74]:


def regression_geneAndExonPositions(genome, k, key, ref = 'MN908947'):
    
    def updateDictionary(Dict, elementList, currentPos, errorList):
        for element in elementList:
            geneOrExon = element.split("|")[0]
            cigar = element.split("|")[2]
            # Empty cigar string
            if cigar == "":
                errorList.append(element)
                continue
            runlen = runLengthEncoding(genome.processCigar(cigar))
            positions = positionWithCigarString(currentPos, cigar)
            # Only an S was seen
            if positions[0] == -1 and positions[1] == -1:
                errorList.append(elementList)
            if geneOrExon in Dict:
                Dict[geneOrExon]["end"] = positions[1]
            else:
                Dict[geneOrExon]["start"] = positions[0] 
                Dict[geneOrExon]["end"] = positions[1]
        return Dict, errorList
    
    # Return the path for the given strain (ABCDEFGH)
    def entityPath(genome, strain=ref):
        path = []
        for i in genome.outgoing["SOURCE"]:
            edge = genome.edges[i]
            # Get the first B edge
            if strain in edge["strain"]:
                current = i 
                path.append(current)
        while current != "SINK":
            # Floating node
            if current[0] == "F" or current[0] == "B":
                for i in genome.outgoing[current]:
                    edge = genome.edges[i]
                    if strain in edge["strain"]:
                        current = i
                        path.append(current)
            # Anchor node
            elif current[0] == "A":
                for edge in genome.outgoing[current]:
                    cEdge = genome.edges[edge]
                    if strain in cEdge["strain"]:
                        current = edge
                        path.append(current)
            # Edge
            else:
                current = genome.edges[current]["dst"]    
                path.append(current)
        # Add the sink node into the path
        if path[-1] == "SINK":
            path = path[:-1]
        for i in genome.incoming["SINK"]:
            edge = genome.edges[i]
            if strain in edge["strain"]:
                path.append(i)
                return path
        return path
    
    def positionWithCigarString(pos, cigar):
        import re
        """
        Given a linear position corresponding to a node and a cigar string representing gene annotations,
        return the position where it starts or stops within the cigar string.
        Cigar should have already been run through the processCigar method and the runLengthEncoding method.
        """
        # Initialize regex patterns
        sms = re.compile(".*S.*M.*S")
        sm = re.compile(".*S.*M")
        ms = re.compile(".*M.*S")
        m = re.compile(".*M")
        s = re.compile(".*S")

        splitString = re.split("S|M", cigar)[:-1]
        for i in range(len(splitString)):
            splitString[i] = int(splitString[i])
        # SMS 
        if(re.search(sms, cigar)):      
            begin = pos + splitString[0] 
#                 begin = pos + splitString[0] - 1
            end = begin + splitString[1] -1
            return (begin, end)
        # SM 
        elif(re.search(sm, cigar)):
            begin = pos + splitString[0]
#                 begin = pos + splitString[0] - 1
            end = begin + splitString[1] - 1
            return (begin, end)
        # MS 
        elif(re.search(ms, cigar)):
            end = pos + splitString[0] - 1
            return (pos, end)
        # M 
        elif(re.search(m, cigar)):
            end = pos + splitString[0] -1
            return (pos, end)
        # S - Error annotation, should report errors here because these edges shouldn't exist
        elif(re.search(s, cigar)):
            return (-1, -1)
              
    # Run length encoding
    def runLengthEncoding(string):
        if len(string) == 0:
            return
        stack = [[string[0], 1]]
        for char in string[1:]:
            if char != stack[-1][0]:
                tupleForChar = [char, 1]
                stack.append(tupleForChar)
            else:
                stack[-1][1] += 1
        output = ""
        for i in stack:
            if i[1] == 1:
                output += i[0]
            else:
                output += (str(i[1])+i[0])
        return output
    
    # Iterate through the genomes checking the gene positions for errors
    def checkGenePositions():
        # Initialize dictionaries and starting position
        nodes = genome.nodes
        edges = genome.edges
        seenGenes = defaultdict(dict)
        
        founder = []
        edgelist = genome.outgoing['SOURCE']
        for edge in edgelist:
            founder += genome.edges[edge]['strain']
        
        currentPos = genome.Nodecoordinates("SOURCE", strainlist = founder)[ref]
        errorList = []
        
        # Get the "B" path from start to end
        path = entityPath(genome, ref)
        for i in range(len(path)):
            # Current node/edge that you are on in the "B" path
            element = path[i]
            # Anchors
            if element[0] == "A":
                # Get the current position so that you set the correct value 
                currentPos = genome.Nodecoordinates(element, strainlist = founder)[ref]
                if key in nodes[element]:
                    genelist = nodes[element][key]
                    updateDictionary(seenGenes, genelist, currentPos, errorList)
            # Floating nodes
            elif element[0] == "F" or element[0] == "B" or element == "SOURCE" or element == "SINK":
                #print element
                currentPos = genome.Nodecoordinates(element, strainlist = founder)[ref]
            # K, L, S, E -> All edge nodes
            else:
                src_node = edges[element]['src']
                if src_node.startswith('A'):
                    currentPos = genome.Nodecoordinates(src_node, strainlist = founder)[ref] + k
                elif src_node == "SOURCE":
                    currentPos = genome.Nodecoordinates(src_node, strainlist = founder)[ref] + 1
                if key in edges[element]:
                    genelist = edges[element][key]
                    updateDictionary(seenGenes, genelist, currentPos, errorList)
        return seenGenes, errorList
            
    sg, el = checkGenePositions()
    return sg, el
    


# In[75]:


keylist = ['CDS', 'five_prime_UTR', 'gene', 'region', 'three_prime_UTR']
k = 11
for key in keylist:
    print( key)
    sg, el = regression_geneAndExonPositions(graph, k, key)
    #print sg
    data = get_gff_data("GCF_009858895.2_ASM985889v3_genomic.gff.gz")
    P1 = parse_feature(data)
    P1
    for i, b in sg.iteritems():
        print( i, b, P1[i]['start'], P1[i]['end'])
        print( b['start'] - P1[i]['start'], b['end'] - P1[i]['end'])


# In[76]:


def get_gff_data(filename):
    with gzip.open(filename, 'rb') as fp:
        data = fp.readlines()
        return data
data = get_gff_data("GCF_009858895.2_ASM985889v3_genomic.gff.gz")
Parser = {}
count = 0
Type = set()
for line in data:
    values = line.split('\t')
    if len(values) < 9:
        continue
    start, end = int(values[3]), int(values[4])
    # parse attributes
    attributes = values[8]
    Att = {}
    attributes_info = attributes.split(';')
    for item in attributes_info:
        k, v = item.split('=')
        Att[k] = v
    #assert "ID" in Att.keys()
    # build in   
    name = Att["ID"]
    Parser[name] = {}
    Parser[name]["start"] = start
    Parser[name]['end'] = end
    Parser[name]['att'] = Att
    Parser[name]['strand'] = values[6]
    Parser[name]['type'] = values[2]
    Type.add(values[2])
Type


# In[77]:


pangenomesize = 0
nodelist = graph.nodes.keys()
edgelist = graph.edges.keys()
pangenomesize += 10 * len(nodelist)
for edge in edgelist:
    pangenomesize += len(graph.edges[edge]['seq'])


# In[78]:


#print pangenomesize, len(sequence)
def contig_name(headers):
    X = []
    for item in headers:
        X.append(item.split(' ')[-1])
    return X

samples = contig_name(header)
#Anchorlist.columns, samples


# In[79]:


total = 0
seqname = list(Anchorlist.columns)
index = [i for i, name in enumerate(samples) if name in seqname]
print len(index)
for i in index:
    seq = sequence[i]
    #print len(seq)
    total += len(seq)
print total, pangenomesize/float(total)


# In[80]:


len(graph.nodes.keys()), len(graph.edges.keys())


# In[81]:


edgenum = numpy.zeros(52)
node = "SOURCE"
while node != "SINK":
    edgelist = graph.outgoing[node]
    edgenum[len(edgelist)] += 1
    node = graph.nextAnchor(node)


# In[82]:


edgenum.astype(int)


# In[83]:


import matplotlib.pyplot as plot
plot.bar(range(25), numpy.log10(edgenum[:25]+1))
plot.xticks(range(25))
plot.xlabel('Haplotype Num')
plot.ylabel("log10(Number of regions)")
plot.show()


# In[84]:


edgenum = numpy.zeros(52)
node = "SOURCE"
while node != "SINK":
    edgelist = graph.outgoing[node]
    if len(edgelist)>10:
        print node
    node = graph.nextAnchor(node)


# In[ ]:




