# SARS-COV-2 Pangenome

We presented a SARS-CoV-2 graphical genome as a reference pangenome, integrating 179 partial and complete sequence (NCBI virus platform) assembly in a single graph representation (updated in time). Here, we generalized the anchor-based linear genome registration and pangenome construction pipeline to construct a SARS-CoV-2 graphical genome. The introduction of anchors, defined as the conserved, unique and consistently ordered sequences in every linear assembly, partitioned and compressed multiple linear genomes in a consistent manner. Sequences between anchors were extracted and merged into edges, which contained all variants in the genome. 


## Genome Registration and Pangenome Construction pipeline applied to SARS-COV-2

We propose both a general methodology for genome registration and pangenome construction based on anchor sequences.
In the previous work, we apply the method to an instance of a mouse genetic reference pangenome, the Collaborative Cross Graphical Genome (CCGG). 
Considering the recent outbreak of the COVID-19, we further apply this method to the SARS-Cov-2 pangenome construction, 
which shows advantage on fast and incrementally integrating accumulated assembly data to a single graph representation. 

## Input
Complete SARS-Cov-2 assemblies were obtained from 
https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049

## Procedure
### 1. Anchor Candidates selection.(AnchorCandidateSelection.py)
The Initial anchor candidates set is population-specific. 
For SARS-Cov-2, currently we select every non-overlapping kmer in the reference assembly as the inital anchor candidates set.
Other options of selecting the initial candidates set so as to reduce the reference bias as much as possible:
1) Go through every existed linear genome and collect all possible kmers from each linear genome, then select the maximum set of the kmer for mapping. 
2) Construct multi-string Burrows–Wheeler transform (msbwt) for every sample reads and find out the maximum set of kmer from sample msbwt.

Noted that options above may lead to a large size of initial mapping candidates. When incrementally adding the linear genome, the new kmers may have to be remapped to the original linear genomes. This may affect the efficiency of genome registration, especially confronting with a fast growing number of collected sequence data.

A potential way to resolve this reasonably is to focus on a set of old samples, i.e. samples collected at the early time. 

Each kmer is of a unique name embedding their relative orders, e.g. the linear coordinates of the reference genome.

### 2. Genome Registration by anchor candidates (genome_registration_2020.py)
We map the anchor candidates back to each individual linear genome. 
We recorded the mapping position for each anchor candidates in a matrix, where the row names is the anchor candidate names 
and column names are the sample name of each linear genome.
If the mapping position is not unique or unordered relative to other anchors in a linear genome, the matrix element will be denoted by "?"

### 3. Dynamic Pangenome Construction Pipeline (genome_registration_2020.py)
Specifying a set of linear genomes that will be incorporated to the pangenome, the mapping position of the anchor candidates 
on these linear genomes can be obtained from the registration table.
The final anchor set can be selected by the anchor candidates that are not "?" in every selected column.
We merge the adjacent anchors by keeping the start and the end of a continuous run.
Each linear genome is further partitioned by the final anchor set. 
Sequences between adjacent anchor pairs were extracted. Identical sequences were merged into an edge with strain annotations.
The graphical pangenome, i.e. anchor node and edge information, is recorded in a fasta file format. Each linear genome can be extracted/recovered from the pangenome through graph traversal.

### 4. Pangenome Annotation(genome_registration_2020.py)
Biological features obtained from the reference gff files were annotated to the reference path of the pangenome.
We further performed the pairwise alignment between each paralleled reference and alternation edge (edges that share the same source and destination node).
Since graph genome compress sequence in a graph structure, the computation load for alignment is significantly reduced.

### 5. Graphical Genome Loading and Variants analysis (CCGG_extension.py, Variants_analysis.py)
The constructed SARS-Cov_2 graphical genome can be loaded and edited by a generlized version of the CCGG API, CCGG_extension.py.
Sequence variants can be obtained by comparing the alignment results and merged into a variants table for futher phologenetic analysis.


