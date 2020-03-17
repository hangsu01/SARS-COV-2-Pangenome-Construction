#!/usr/bin/env python
# coding: utf-8

# In[2]:


import CCGG_extension as CCGG
import pandas as pd
import time
import numpy
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


# In[4]:


start = time.time()
chromo = '1'
#graph = CCGG.GraphicalGenome("Nodefile.fa","Edgefile.fa", nnamelen= 6, enamelen=7)
graph = CCGG.GraphicalGenome("CNCB_Cov_Nodefile.fa","CNCB_Cov_Edgefile.fa", nnamelen= 6, enamelen=7)

print "Reading graph %4.2f secs" % (time.time() - start)


# In[1]:


def get_variant_position(cigar):
    ref_pos = []
    alt_pos = []
    alt_i = 0
    ref_i = 0
    for i, s in enumerate(cigar):
        if s == 'I':
            if ref_i > 0:
                ref_pos.append(ref_i-1)
            else:
                ref_pos.append(ref_i)    
            alt_pos.append(alt_i)
            alt_i += 1
        if s == 'D':
            ref_pos.append(ref_i)
            if alt_i > 0:
                alt_pos.append(alt_i-1)
            else:
                alt_pos.append(alt_i)
            ref_i += 1
            
        if s == 'X':
            ref_pos.append(ref_i)
            alt_pos.append(alt_i)
            alt_i += 1
            ref_i += 1
            
        if s == '=':
            alt_i += 1
            ref_i += 1

    return ref_pos, alt_pos


# In[7]:


founder = []
for edge in graph.outgoing['SOURCE']:
    founder += graph.edges[edge]['strain']
print len(founder)


# In[8]:




refstrain = 'MN908947'
k = 10

def find_allVar(graph, refstrain, founder, k):
    Var = {} # Var[refpos]['edgename'] = ['A']
    node = "SOURCE"
    while node != "SINK":
        edgelist = graph.outgoing[node]
        
        for edge in edgelist:
            if refstrain in graph.edges[edge]['strain']:
                refseq = graph.edges[edge]['seq']
                refedge = edge
                break
                
        for edge in edgelist:
            cigar = graph.edges[edge]['variants']
            if node != "SOURCE":
                refstart = int(graph.Nodecoordinates(node, strainlist = founder)[refstrain]) + k
            else:
                refstart = int(graph.Nodecoordinates(node, strainlist = founder)[refstrain])
            
            refpos, altpos = get_variant_position(graph.processCigar(cigar))
            #print refpos, altpos
            alt_seq = graph.edges[edge]['seq']
            for i, rp in enumerate(refpos):
                pos = refstart + rp
                Var[pos] = Var.get(pos, {})
                Var[pos][edge] = Var[pos].get(edge, {})                
                Var[pos][edge]['base'] = Var[pos][edge].get('base', "") + alt_seq[altpos[i]]
                Var[pos][edge]['altoffset'] = Var[pos][edge].get('altoffset', []) + [altpos[i]]
                Var[pos][refedge] = Var[pos].get(refedge, {})
                Var[pos][refedge]['base'] = refseq[refpos[i]]
        node = graph.nextAnchor(node)    #print cigar, refpos, altpos
    return Var


# In[9]:


def coordinate_exchange(Graph, coord, alter_strain, refstrain, k = 10):
    def find_offset_on_path(Graph, itemlist, entity, offset):
        prevoff = 0
        for item in itemlist:
            if item == entity:
                path_coord = prevoff + offset
                return path_coord
            else:
                if item.startswith('A'):
                    prevoff += k
                elif item.startswith('F') or item.startswith('B') or item == 'SOURCE' or item == 'SINK':
                    prevoff += 0
                else:
                    prevoff += len(Graph.edges[item]['seq'])

    def completecigar(Graph, edgelist, kmer = k):
        cigar = ''
        for edge in edgelist:
            if edge.startswith('B') or edge.startswith('F') or edge == 'SOURCE' or edge == 'SINK':
                continue
            if edge.startswith('A'):
                cigar += str(kmer) + '='
                continue
            if 'variants' in Graph.edges[edge].keys():
                cigar += Graph.edges[edge]['variants']
            else:
                #print edge
                return False
        else:
            return cigar

    def findsequence(Graph, pathlist):
        seq = ''
        for item in pathlist:
            if item.startswith('A'):
                seq += '' # do not count anchor length
            elif item.startswith('L') or item.startswith('E')or item.startswith('K'):
                seq += Graph.edges[item]['seq']
            elif item.startswith('S') and item != "SOURCE" and item != 'SINK':
                seq += Graph.edges[item]['seq']
            else:
                seq += ''
        return seq
    
    def offsetexchange(cigar, B_offset):
        alt_i = 0
        ref_i = 0
        for i, s in enumerate(cigar):
            if ref_i == B_offset:
                return alt_i   
            if s == '=':
                alt_i += 1
                ref_i += 1

            if s == 'I':
                alt_i += 1
            if s == 'D':
                ref_i += 1
            if s == 'X':
                alt_i += 1
                ref_i += 1
            
    def Complete_Parallel_edge_offset_exchange(Graph, coord, alter_strain, founderlist, refstrain, output = 0):
        """Alignment-based entity offset exchange
        (only apply to parallel edge, edge share the same src and dst, offset exchange in Path scale are not finished yet)
        Given entity+offset on B path, return the position on the alternative strain
        If alignment results are not applicable, return None

        Parameters:
            entity: <str> - Graph entity on B path
            offset: <int> - offset on the entity
            alter_strain: <str> - strain attributes for the target alternative path position
            strain: <default> "B" path, when multi-alignment cigar are added, this could be further implemented
        """
        entity, offset = Graph.linear_to_entityoffset(coord, refstrain)

        if entity.startswith('A'):
            alter_coord = Graph.Nodecoordinates(entity, strainlist=founderlist)[alter_strain] + offset
            return alter_coord, offset
        else:
            src = Graph.edges[entity]['src']    
            if src.startswith('A') or src == 'SOURCE':
                s_anchor = src
            else:
                s_anchor = Graph.prevAnchor(src)
            d_anchor = Graph.nextAnchor(s_anchor)

            # construct path
            Path = Graph.getSubPath(s_anchor, d_anchor, init_strain=founderlist)
            for itemlist, strain in Path:
                if refstrain in strain and alter_strain in strain:
                    path_coord = find_offset_on_path(Graph, itemlist, entity, offset)
                    alter_coord = Graph.Nodecoordinates(s_anchor, strainlist = founderlist)[alter_strain] + path_coord
                    return alter_coord, path_coord
                
                elif refstrain in strain:
                    B_offset = find_offset_on_path(Graph, itemlist, entity, offset)
                    
                elif alter_strain in strain:
                    c_cigar = completecigar(Graph, itemlist)
                    if c_cigar == False: # not applicable in covid19 pangenome
                        for il, bstrain in Path:
                            if refstrain in bstrain:
                                ref_seq = Graph.findsequence(il, countinganchor=True)
                                break
                        alt_seq = Graph.findsequence(itemlist, countinganchor=True)
                        cigar = makeCigar(alt_seq, ref_seq)
                        cigar = Graph.processCigar(cigar)
                    else:
                        cigar = Graph.processCigar(c_cigar)

            path_coord = offsetexchange(cigar, B_offset)
            #print cigar, len(cigar),B_offset
            alter_coord = Graph.Nodecoordinates(s_anchor, strainlist=founderlist)[alter_strain] + path_coord
            return alter_coord, path_coord
        
    founderlist = []
    for edge in Graph.outgoing["SOURCE"]:
        founderlist += Graph.edges[edge]['strain']
    
    alt_coord, path_coord = Complete_Parallel_edge_offset_exchange(Graph, coord, alter_strain, set(founderlist), refstrain)
    
    return alt_coord, path_coord


# In[11]:


Var = find_allVar(graph,refstrain, founder, k)
#Var


# In[12]:


# collapse adjacent 
def collapsed_adjacent(Hap_Dict, step):
    k = step
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
        pos_inblocks = []
        for s,e in blocks:
            pos_inblocks.append((poslist[s], poslist[e-1]))
            
        for s,e in blocks:
            for i in range(s,e):
                errorindex.append(i)
                
        return errorindex, pos_inblocks

    RefIndexlist = sorted(Hap_Dict.keys())
    errorindex, pos_inblocks = sort_by_kmerlength(RefIndexlist, k)
    index = [RefIndexlist[i] for i in range(len(RefIndexlist)) if i not in errorindex]
    index += pos_inblocks
    
    return index


# In[3]:


def Var2Hap(Var, founder):
    RefIndex = collapsed_adjacent(Var, 1)
    
    Hap = {} # Hap[refpos][haplotype] = [edges]
    Ns = set()
    #print len(RefIndex)
    for refpos in RefIndex:
        if isinstance(refpos, int):
            Var_D = Var[refpos]
        else:
            s,e = refpos # both end included
            Var_D = {}
            sanchor, eanchor = graph.boundingAnchors(refstrain, s,e)
            Path = graph.getSubPath(sanchor, eanchor, init_strain=set(founder))
            
            for itemlist, strains in Path:
                edge = itemlist[1]
                assert edge.startswith('E') or edge.startswith("S")

                sequence = graph.findsequence(itemlist,countinganchor=True)
                alter_strain = list(strains)[0]
                alts, soffset = coordinate_exchange(graph, s, alter_strain, refstrain)
                alte, eoffset = coordinate_exchange(graph, e, alter_strain, refstrain)
                haplotype = sequence[soffset:eoffset+1]
                Var_D[edge] = {}
                Var_D[edge]['altoffset'] = [range(soffset,eoffset+1)]
                Var_D[edge]['base'] = haplotype
        H = {}
        for edge, D in Var_D.iteritems():
            haplotype = D['base']
            
            if "N" in haplotype: # exlude all Ns
                Ns.add(edge)
                continue
            
            altpos = D.get('altoffset', "?")
            if altpos == "?":
                H[haplotype] = H.get(haplotype, []) + [(edge, refpos)]
            else:
                H[haplotype] = H.get(haplotype, []) + [(edge, altpos)]
            
        # check 
        if len(H)>1:
            Hap[refpos] = H
        else:
            print H, refpos
    return Hap, Ns

Hap, Nedge = Var2Hap(Var, founder)
# print len(Hap)
# print Nedge
#Hap[(21627, 21707)]
#print Var[0]
#Var2Hap(Var)


# In[14]:


#print Hap[Hap.keys()[0]]

# create Haplotype table, position by sample, entry is the sequence 
def get_haplotype_table(graph, Hap):
    Info = {}
    for refpos, Var_D in Hap.iteritems():
        for haplotype, items in Var_D.iteritems():
            for edge, altpos in items:
                strainlist = graph.edges[edge]['strain']
                for strain in strainlist:
                    Info[refpos] = Info.get(refpos, {})
                    Info[refpos][strain] = haplotype
    df = pd.DataFrame(Info)
    df = df.dropna(axis=1) # delete position with ambiguous sequences
    return df
df = get_haplotype_table(graph, Hap)
print df.shape
df.to_csv('CovVariants.csv')

# # # drop the prefix and suffix of the graphical genome
# df = df.drop([(0,33), (29827, 29903)], axis=1)
# print df.shape


# In[15]:


# construct levenshtein distance matrix
from Levenshtein import distance

def construct_distancemat(df, excludesample = []):
    excludesample = [unicode(item) for item in excludesample]
    df = df.drop(excludesample, axis=0)
    s, p = df.shape
    print s, p
    labels = df.index
    dist = numpy.zeros((s,s))
    #print dist
    for i in range(s):
        h0 = df.iloc[i, :].values
        for j in range(s):
            h1 = df.iloc[j, :].values
            d = 0
            for t in range(p):
                #print h0[t], h1[t]
                if h0[t] == h1[t]:
                    continue
#                 if len(h0[t]) == len(h1[t]):
#                     d += sum([1 for x, y in zip(h0[t], h0[t]) if x.upper() != y.upper()])
#                 else:
                d += distance(h0[t], h1[t])
            dist[i,j] = d
    return dist, labels
dist, labels = construct_distancemat(df, excludesample= ['NMDC60013002-05'])
print dist.shape


# In[16]:


#df.drop(u'NMDC60013002-05', axis = 0)


# In[18]:


# get sample info
table = 'Cov_Sample_info.csv'
sampleinfo = pd.read_excel(table, header=1, index_col=None)
sampleinfo
#print sampleinfo.head()
def isNaN(num):
    return num != num

SampleD = {}
for i in range(len(sampleinfo)):
    sample = sampleinfo['Accession ID'][i]
    region = sampleinfo['Location'][i]
#     if isNaN(region):
#         region = sampleinfo.iloc[i,:]['GenBank_Title'].split(',')[0]
    SampleD[sample] = region


# In[19]:


print labels
l = [SampleD[str(i)] for i in labels]
print l
#SampleD


# In[20]:


D = pd.DataFrame(dist, columns = labels, index = l)
D.to_csv('CovDistancematrix.csv')


# In[21]:


from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
Z = hierarchy.linkage(data_scaled, method="average") # method = "average" UPGMA; "weighted", WPGMA; ’centroid’ WPGMC
plt.figure(figsize=[20,10])
l = [SampleD[i] for i in labels]
#l = [i.split('/')[0] for i in l]
print l
dn = hierarchy.dendrogram(Z,labels=labels,leaf_font_size=13)
dn['ivl'][0]


# In[22]:


sample = 'NMDC60013002-05'
node = "SOURCE"
counter = numpy.zeros(100)
while node != 'SINK':
    edgelist = graph.outgoing[node]
    for edge in edgelist:
        if sample in graph.edges[edge]['strain']:
            counter[len(graph.edges[edge]['strain'])] += 1
            break
    node = graph.nextAnchor(node)
print counter, sum(counter), len(graph.nodes.keys()), 84/381.0


# In[23]:


# filter position without variants (according sample set)
#print "df", df.shape


# In[24]:


mydata[:,100]


# In[25]:


df.index


# In[191]:


from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# colors = 'lime red blue magenta yellow'.split()
# cmap = matplotlib.colors.ListedColormap(colors, name='colors', N=None)
# print cmap
# plt.imshow(p, cmap=cmap)
# plt.show()
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable


def Plot_Indels(df,roworder = df.index, shift = -0.5, excludesample = []):
    '''Import the pandas dataframe where columns correspond to ref intervals, rows correspond to strains, 
    elements denote haplotypes
     roworder is the order for samples, to get a pretty plot, input the clustered order
     
    '''
    excludesample = [unicode(item) for item in excludesample]
    roworder = [item for item in roworder if item not in excludesample]
    def filternonvariantsposition(df, roworder):
        data = df.loc[roworder]
        colorder = data.columns

        # filter position without variants (according sample set)
        delcolumns = []
        for i in range(len(colorder)):
            if len(set(data.iloc[:,i].values)) < 2:
                delcolumns.append(colorder[i])
        colorder = [item for item in colorder if item not in delcolumns]


        data = data.loc[:,colorder]
        return data

    data = filternonvariantsposition(df, roworder)
    print "data", data.shape

    def expandforplot(data):
        columns = []
        for i in data.columns:
            if isinstance(i, int):
                columns.append(i)
            else:
                columns.append(i[0])

        index = numpy.argsort(columns)
        colorder = [data.columns[i] for i in index]
        data = data.loc[:,colorder]
        p = data.values            


        def expand_indels(index, p):
            r,c = p.shape
            #print p[:,index]
            m = [len(seq) for seq in p[:,index]]
            m = max(m)
            a = numpy.ndarray([r,m], dtype=object)

            for i, seq in enumerate(p[:,index]):
                for j in range(len(seq)):
                    a[i,m - len(seq) + j] = seq[j]
            return a



        def expandmatrix(p):
            mydata = expand_indels(0, p)
            labels_y = [0] + [""] * (mydata.shape[1]-2) + [33]
            r, c = p.shape
            for index in range(1,c):
                index_ = colorder[index]
                if isinstance(index_, int):
                    mydata = numpy.hstack((mydata, p[:,index].reshape(-1,1)))
                    assert len(set(p[:,index]))>1
                    labels_y += [index_]
                else:
                    a = expand_indels(index, p)
                    mydata = numpy.hstack((mydata, a))
                    labels_y += [index_[0]] + [""] * (a.shape[1]-2) + [index_[1]]
            return mydata, labels_y

        mydata, labels_y = expandmatrix(p)

        return mydata, labels_y

    mydata, labels_y = expandforplot(data)
    print mydata.shape, labels_y
    
    
    # exchange character to number for plot
    character = set()
    for i in mydata:
        for j in i:
            character.add(j)

    colorD = {}
    count = 1
    for i in ['A', 'C', 'G', 'T', 'S', 'R', 'W', 'Y']:
        colorD[i] = count
        count += 1

    r,c = mydata.shape
    p = mydata
    print p.shape
    
    for i in range(r):
        for j in range(c):
            p[i,j] = colorD.get(p[i,j], 0)
    p = p.astype(float)
    
    # plot
    colorDict = {'A': 'blue', "G":'green', 'C':'red', 'T': 'yellow', "None": 'grey'}
    cmap = mpl.colors.ListedColormap(["white",'blue', 'red', 'green', 'yellow', 'cyan', 'pink', 'grey','orange'], name = 'colors', N = None)
    n = mpl.colors.Normalize(vmin=0,vmax=3)

    plt.figure(figsize=[200,400])
    AX = plt.gca()
    im = AX.matshow(p, cmap=cmap)
    plt.xlim([shift-0.5,c])
    plt.ylim([shift-0.5,r])
    # set ticks
    plt.yticks(range(0,r), roworder,rotation=0, fontsize = 10)
    plt.xticks(range(0,c), labels_y, rotation = 80, fontsize = 10)

#     #plot grid
    r,c = p.shape
    
    ax = np.arange(shift,c,1)
    for z in ax:
        plt.plot((z,z),(shift ,r+shift),'k-')

    ax = np.arange(shift,r,1)
    for z in ax:
        plt.plot((shift,c+shift),(z,z),'k-')

#     #display bases
#     data = df.loc[roworder]
#     data = data.loc[:,colorder]
#     Bases = data.values
    
#     for i in range(r):
#         for j in range(c):
#             base = Bases[i,j]
#             AX.text(j, i, str(base), va='center', ha='center')

    #legend
    #cbar = plt.colorbar(im)
    character = ["",'A', 'C', 'G', 'T', 'S', 'R', 'W', 'Y']
    divider = make_axes_locatable(AX)
    cax = divider.append_axes("right", size="1%", pad=0.005)

    cbar = plt.colorbar(im, cax=cax)

    cbar.ax.get_yaxis().set_ticks([])

    # annotate characters
    for j, lab in enumerate(character):
        cbar.ax.text(.5, (j+0.5)/9.0, lab, ha='center', va='center',fontsize = 10)
    # cbar.ax.get_yaxis().labelpad = 15
    # cbar.ax.set_ylabel('# of contacts', rotation=270)

    # save
    plt.show()
    plt.savefig("Variantsdisply.png")


# In[190]:


roworder = dn['ivl']
print roworder
Plot_Indels(df,roworder)


# In[27]:


from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_snps(df,roworder = df.index, shift = -0.5, excludesample = []):
    excludesample = [unicode(item) for item in excludesample]
    roworder = [item for item in roworder if item not in excludesample]
    
    def filternonvariantsposition(df, roworder):
        data = df.loc[roworder]
        colorder = data.columns

        # filter position without variants (according sample set)
        delcolumns = []
        for i in range(len(colorder)):
            if len(set(data.iloc[:,i].values)) < 2:
                delcolumns.append(colorder[i])
        colorder = [item for item in colorder if item not in delcolumns]


        data = data.loc[:,colorder]
        return data

    data = filternonvariantsposition(df, roworder)
    print "data", data.shape
    
    def sortcolumns(data):
        columns = []
        for i in data.columns:
            if isinstance(i, int):
                columns.append(i)
        data = data.loc[:,columns]        
        return data, columns
    
    data,columns = sortcolumns(data)
    p = data.values
    
    character = set()
    for i in data.values:
        for j in i:
            character.add(j)
    print character

    colorD = {}
    count = 1
    for i in ['A', 'C', 'G', 'T', 'N']:
        colorD[i] = count
        count += 1

    r,c = data.shape
    for i in range(r):
        for j in range(c):
            if p[i,j] not in "ACGT":
                p[i,j] = colorD['N']
            else:
                p[i,j] = colorD[p[i,j]]
    p = p.astype(float)


    colorDict = {'A': 'blue', "G":'green', 'C':'red', 'T': 'yellow'}
    cmap = mpl.colors.ListedColormap(['blue', 'red', 'green', 'yellow', "grey"], name = 'colors', N = None)
    #n = mpl.colors.Normalize(vmin=0,vmax=3)

    plt.figure(figsize=[20,20])
    AX = plt.gca()
    im = AX.matshow(p, cmap=cmap)
    plt.xlim([-0.5,c-0.5])
    plt.ylim([-0.5,r-0.5])
    # set ticks
    # print roworder
    # print columns
    plt.yticks(range(0,r), roworder,rotation=0)
    plt.xticks(range(0,c), columns, rotation = 50)

    #plot grid
    r,c = p.shape
    ax = np.arange(-0.5,c,1)
    for z in ax:
        plt.plot((z,z),(-0.5,r-0.5),'k-')

    ax = np.arange(-0.5,r,1)
    for z in ax:
        plt.plot((-0.5,c-0.5),(z,z),'k-')

    # display bases
    data = df.loc[roworder]
    data = data.loc[:,columns]

    for i in range(r):
        for j in range(c):
            base = data.iloc[i,j]
            AX.text(j, i, str(base), va='center', ha='center')

    #legend
    #cbar = plt.colorbar(im)
    character = ['A', 'C', 'G', 'T', 'N']
    divider = make_axes_locatable(AX)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    cbar = plt.colorbar(im, cax=cax)

    cbar.ax.get_yaxis().set_ticks([])


    for j, lab in enumerate(character):
        cbar.ax.text(.5, (j+0.5)/5.0, lab, ha='center', va='center',fontsize = 20)
    # cbar.ax.get_yaxis().labelpad = 15
    # cbar.ax.set_ylabel('# of contacts', rotation=270)


    plt.show()
    plt.imsave("SNPSdisply.png", im)

plot_snps(df, roworder = dn['ivl'])


# In[196]:


u'MT019530' in df.index


# In[32]:


# create variants table
roworder = dn['ivl']
columns = []
for i in df.columns:
    if isinstance(i, int):
        columns.append(i)
    else:
        columns.append(i[0])
        
index = numpy.argsort(columns)
colorder = [df.columns[i] for i in index]
#print colorder
df = df.loc[roworder]
df = df.loc[:,colorder]
#df

Hap = {}
for i in range(len(colorder))[1:-1]:
    H = {}
    sequences = df.iloc[:,i].values
    for j in range(len(sequences)):
        s = sequences[j]
        H[s] = H.get(s, []) + [roworder[j]]
    Hap[colorder[i]] = H
for pos, item in Hap.iteritems():
    print pos, item.keys(), [len(item[k]) for k in item.keys()]


# In[377]:


def offsetexchange(cigar, B_offset):
    alt_i = 0
    ref_i = 0
    for i, s in enumerate(cigar):
        #print ref_i
        if ref_i == B_offset:
            return alt_i   
        if s == '=':
            alt_i += 1
            ref_i += 1

        if s == 'I':
            alt_i += 1
        if s == 'D':
            ref_i += 1
        if s == 'X':
            alt_i += 1
            ref_i += 1

def get_variant_position(cigar):
    ref_pos = []
    alt_pos = []
    alt_i = 0
    ref_i = 0
    for i, s in enumerate(cigar):
        if s == 'I':
            if ref_i > 0:
                ref_pos.append(ref_i-1)
            else:
                ref_pos.append(ref_i)    
            alt_pos.append(alt_i)
            alt_i += 1
        if s == 'D':
            ref_pos.append(ref_i)
            if alt_i > 0:
                alt_pos.append(alt_i-1)
            else:
                alt_pos.append(alt_i)
            ref_i += 1
            
        if s == 'X':
            ref_pos.append(ref_i)
            alt_pos.append(alt_i)
            alt_i += 1
            ref_i += 1
            
        if s == '=':
            alt_i += 1
            ref_i += 1

    return ref_pos, alt_pos

print graph.edges['S000001']['seq']           
cigar = graph.processCigar('2I1X3=2I2=2X3=1D2=3X23=')
print offsetexchange(cigar, 0)

print get_variant_position(cigar)


# In[693]:


#plot
y,x = numpy.histogram(hapnum, bins = len(conseq)/500, range = [0, len(conseq)])
# print x
# print y
plt.plot(x[:-1], y, color='b', linestyle='-', alpha=0.5, linewidth = 1.5)
plt.xlabel('Genome/kb')
plt.ylabel('Haplotype Num')


# In[694]:


# delete the prefix and suffix
plt.plot(x[8:-4], y[8:-3])
plt.xlabel('Genome/kb')
plt.ylabel('Haplotype Num')


# In[309]:


# overlay gene
import gzip
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


# In[310]:


data=get_gff_data("GCF_009858895.2_ASM985889v3_genomic.gff.gz")
Parser = parse_feature(data)


# In[311]:


gene_d = {}
interval = []
for name,D in Parser.iteritems():
    if D['type'] == 'gene':
        gene_d[name] = {}
        gene_d[name]['start'] = D['start']
        gene_d[name]['end'] = D['end']
        interval.append((D['start'], D['end']))
gene_d


# In[312]:


colors = ['r','g','b','k','grey','orange','purple','pink','darkgreen','darkblue','darkred']
plt.plot(x[2:-1], y[1:-1])
#plt.axhline(y=20, xmin=0, xmax= 1, color = 'red', linestyle='-', linewidth=5)
i = 0
interval = sorted(interval)
for s, e in interval:
    print s, e,
    print float(s)/(3*10**4), float(e)/(3*10**4)
    plt.axhline(y=max(y[1:-1])+1+(i/2.0), xmin=float(s)/(3*10**4), xmax= float(e)/(3*10**4), 
                color = colors[i], linestyle='-', linewidth=2)
    i = i + 1
plt.xlabel('Genome/kb')
plt.ylabel('Haplotype Num')
plt.legend(["haplotype distribution","gene"], loc = "best")

    


# In[ ]:





# In[ ]:





# In[ ]:




