# SARS-COV-2 Pangenome
The emergence of infectious disease, COVID-19, caused by the severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2, previously named as 2019-nCoV) has caused a pandemic and becomes a major threat to global public health. As of March 23, data from John Hopkins Coronavirus resource center (Dong et al. 2020) has shown that more than 350,000 cases have been confirmed in 167 countries and regions. Despite the mild symptoms in most cases, disease onset may result in progressive respiratory failure and cause deaths. COVID-19 is an emerging, rapid evolving situation. The future evolution, adaption and spread of the SARS-CoV-2 warrant urgent monitoring. SARS-CoV-2 belongs to the Coronaviridae family, possessing a single-strand, positive-sense RNA genome, ranging from 26-32 kilobases (Lu et al. 2020). It is the seventh coronavirus know to infect humans (Corman et al. 2018). Whole-genome characterization suggests the SARS-CoV-2 genome shares 79.6% sequence identity with SARS-CoV and is 96% identical to a bat coronavirus, RaTG13 (Wu et al. 2020). The amino acid sequences of the seven conserved replicase domain in ORF1ab, which served as the identity of CoV species, were 94.4% identical between SARS-CoV2 and SARS-CoV, suggesting that the two viruses belong to the same species (Zhou et al. 2020). Due to high error rates of RNA replications, RNA viruses exist as quasi-species, presenting as an ensemble of genotypes with large numbers of variants (Bull et al. 2005). Genetic diversity within a quasi-species has been proposed to contribute to pathogenesis (Denison et al. 2011). Investigating the quasi-species dynamics and viral fitness benefits for deciphering the molecular basis for phenotypic flexibility, with implications for antiviral therapy.
The reference genome of SARS-CoV-2 is desgnated as NC_045512.2 in NCBI Reference Sequence and GeneBank accession number MN908947. This linear reference genome is the whole genome sequence (29,903 bp) of the bronchoalveolar lavage fluid sample collected from a single patient who worked at Huanan Seafood Market at Wuhan (Wu et al. 2020). Reference genome served as the basis for downstream SARS-CoV-2 genome characterization and comparison. However, a linear reference genome fails to capture the genetic diversity within a quasi-species population. Due to the limited RNA replication fidelity, the deviation from a genomic sequence collected from a single patient to quasi-species population is not trivial and will introduce reference bias to the downstream analysis. Recently, the concept of pangenome has been proposed and be regarded as the path forward to comprehensively and efficiently represent the genomic diversity within a population. Pangenome is a collection of genomic sequences to be analyzed jointly as a reference (Computational Pan-Genomics 2018). Graph is a commonly used representation for pangenomes. It provides a natural framework for representing shared sequences and variations among them. A comprehensive reference pangenome is also compatible to downstream tools and facilitate computational efficiency.

In this study, we presented a SARS-CoV-2 graphical genome as a reference pangenome, integrating 179 partial and complete sequence assembly in a single graph representation (updated in time). In the previous work, (Su et al. 2019) proposed a pangenome instance, Collaborative Cross Graphical Genome (CCGG), in the Collaborative Cross population. Here, we generalized the method to construct SARS-CoV-2 graphical genome. The introduction of anchors, defined as the conserved, unique and consistently ordered sequences in every linear assembly, partitioned and compressed multiple linear genomes in a consistent manner. Sequences between anchors were extracted and merged into edges, which contained all variants in the genome. 


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


