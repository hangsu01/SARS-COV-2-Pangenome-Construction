# SARS-COV-2 Pangenome
## Genome Registration and Pangenome Construction pipeline applied to SARS-COV-2

The rapid growth of sequencing data across multiple individuals and populations has driven a paradigm shift from linear reference genome assemblies to 
pangenome representations that incorporate multiple genomic sequences.
Reference genome assemblies play an essential role by providing a coordinate framework for referring to and annotating biological and sequence features. 
However, there are few current standards for assigning coordinates for a graph-based pangenome representations with varying loci.

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
For SARS-Cov-2, we select every non-overlapping kmer in the reference assembly as the inital anchor candidates set.
Each kmer is of a unique name embedded by their linear coordinates of the reference genome.

### 2. Genome Registration by anchor candidates (genome_registration_2020.py)
We map the anchor candidates back to each individual linear genome. 
We recorded the mapping position for each anchor candidates in a matrix, where the row names is the anchor candidate names 
and column names are the sample name of each linear genome.
if the mapping position is not unique or unordered relative to other anchors in a linear genome, the matrix element will be denoted by "?"

### 3. Dynamic Pangenome Construction Pipeline (genome_registration_2020.py)
Specifying a set of linear genomes that will be incorporated to the pangenome, the mapping position of the anchor candidates 
on these linear genomes can be obtained from the registration table.
The final anchor set can be selected by the anchor candidates that are not "?" in every selected column.
We merge the adjacent anchors by keeping the start and the end of a continuous run.
The anchors and their mapping position is recorded in the Node fasta files.
Each linear genome is further partitioned by the anchors. 
Sequences between adjacent anchor pairs were extracted. Identical sequences were merged into an edge with strain annotations.
The edge sequences, their source, destination node and strain information is recorded in the edge fasta file.

### 4. Pangenome Annotation(genome_registration_2020.py)
After creating the SARS-Cov-2 pangenome. 
Biological features obtained from the reference gff files were annotated to the reference path of the pangenome.
We further performed the pairwise alignment between each paralleled reference and alternation edge, i.e. edges that share the same source and destination node.
Since graph genome compress sequence in a graph structure, the computation load for alignment is greately reduced.

### 5. Variants analysis ()
Sequence variants can be obtained by comparing the alignment results and merged into a variants table for futher phologenetic analysis.

