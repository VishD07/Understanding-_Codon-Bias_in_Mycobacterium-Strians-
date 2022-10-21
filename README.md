Understanding Codon Bias and Gene Clustering in Mycobacterium species

In this project I took two mycobacterium species one pathogenic ( Mycobacterium tuberculosis) and other non-pathogenic ( Mycobacterium Smegmatis ) 
The project requires CDS( coding sequence region) fasta file.
It usages RSCU values and clustering techniques to identify codon preferences and explores GC content at various codon position. the goal is to understand difference between pathogenic and non pathogenic species.
. Understand the molecular mechanisms and evolutionary basis underlying differences between pathogenic and non-pathogenic species.
. Examine the role of mutation and selection in shaping codon usage patterns.
. Provide insights into the genomic adaptations linked to pathogenicity.

How code works:
Firstly code takes CDS fasta file calculate codon frequency RSCU values then plot neutrality and measures k-value by elbow method and do the clustering


To run the code you must have cds files, all libraries and python.

./python code.py --input path_to_input_file --output path_to_output_file

all charts will generate in output directory

The project creates four outputs 

1. Codon bias and rscu value bar chart by amino acids.
2. Neutrality Plot - A regression graph showing GC3 vs. GC12
3. Elbow curve chart - to determine value of k which will be further applied in clustering 
4. Clustering Graph - Genes grouped on codon bais using PCA and K-means

I have also added requirement file of all libraries which are needed. 
