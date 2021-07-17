# Biological Sequences of SARS-Related Corona Viruses: Data Analysis
## Calculating the **mean** and **standard deviation** of the coronaviruses genomes' lengths and GC content
> working with genomic data of MERS-CoV, SARS-CoV, and SARS-CoV2 collected in Asia, South America, and Africa 
#
> setting up
```julia
using BioinformaticsBISC195

# first, parse the data into 2 vectors, where one contains parsed headers and the other contains entire sequences
genomes = parse_fasta("data/SARS-related-coronaviruses-genomes.fasta") 
```
> length
```julia 
# now that the data is parsed and ready to be analyzed, the mean and standard deviation of the lengths of the coronaviruses' genomes can be determined using some functions from the Statistics package

using Statistics

#mean of lengths
mean_length = mean([length(seq) for seq in genomes[2]])

#standard deviation of lengths
std_length = std([length(seq) for seq in genomes[2]]) 
```
> GC content

```julia
# using the same parsed data and the Statistics package to determine the mean and standard deviation of the coronaviruses' genomes' GC content (excluding ambiguous bases)

#mean of GC content
mean_GC = mean([gc_content(seq) for seq in genomes[2]])

#standard deviation of GC content
std_GC = std([gc_content(seq) for seq in genomes[2]])
```
