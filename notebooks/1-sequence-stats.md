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
> Length: data
```julia 
# now that the data is parsed and ready to be analyzed, the mean and standard deviation of the lengths of the coronaviruses' genomes can be determined using some functions from the Statistics package

using Statistics
# mean of lengths
mean_length = mean([length(seq) for seq in genomes[2]]) #21863.314634698032
# standard deviation of lengths
std_length = std([length(seq) for seq in genomes[2]]) #12728.185835726168
```


> GC content: data

```julia
# using the same parsed data and the Statistics package to determine the mean and standard deviation of the coronaviruses' genomes' GC content (excluding ambiguous bases)

# mean of GC content
mean_GC = mean([gc_content(seq) for seq in genomes[2]]) #0.3872921217727834

# standard deviation of GC content
std_GC = std([gc_content(seq) for seq in genomes[2]]) #0.022013502530253505
```
> Initial analysis

While the values of the mean length and GC content of the coranaviruses' genomes serve an informative purpose, the standard deviations can give us an idea of how related the genomes are.\
A large standard deviation indicates large variances from sample to sample, and a small standard deviation indicates the opposite.

The standard deviation of genome length is **12728.185835726168 basepairs**, which indicates that the genome length differs rather substantially from genome to genome.\
This makes sense considering that the data contains genomes from three different kinds of coronaviruses.

However, the standard deviation of genome GC content is **0.022013502530253505 basepairs**, which indicates that the genomes contain a considerably similar proportion of G and C basepairs.   

It can be crudely concluded that the coronaviruses' genomes are considerably similar and therefore closely related in terms of base composition (GC content) but also vary drastically in lengths. 

> Cleaning the data

The large standard deviation of genome lengths suggests that there are outliers in the dataset.
```julia
using Plots

# making an Array for the genome lengths
genome_lengths = [length(seq) for seq in genomes[2]]

# making a histogram to display counts of genome lengths: 
histogram(genome_lengths, legend = false, xlabel="genome length (# of basepairs)", ylabel="count", nbins=10^5)
```
The cluster of bars on the far left of the graph confirms the prediction above; there are quite a few genomes with lengths less than (approximately) 5000 basepairs.

These outliers must be removed to not skew the data from having a normal distribution.

To achieve this, I will write a new function ```removeUnder25k(parsed_data)``` in the BioinformaticsBISC195 package to take out any genomes shorter than 25k basepairs and display this "cleaned" data with a new histogram.

```julia

# storing the cleaned data in the variable `genomes_cleaned`
genomes_cleaned = removeUnder25k(genomes)

# making an Array for the genome lengths from the cleaned data
genomes_lengths_cleaned = [length(seq) for seq in genomes_cleaned[2]]

# making a new histogram with the cleaned data to display counts of genome lengths:
histogram(genomes_lengths_cleaned, legend = false, xlabel="genome length (# of basepairs)", ylabel="count", nbins=400)
```
Now, the cluster of bars on the far left of the histogram are gone. The distribution of counts of genome lengths is now fairly normal. 