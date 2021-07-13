# BISC 195 Final Project Analysis Proposals 
# Emily Lu

## Preface
Below you will find two of my proposed analysis plans. They will generally follow the following organization and order:

- Main question of analysis
- Analysis ideation (if applicable)
- Sub-questions of analysis

    1. Subquestion of analysis 1
        - A list of Julia function(s) utilized
        - Method(s) for data visualization
    2. Subquestion of analysis 2
        - A list of Julia function(s) utilized
        - Method(s) for data visualization
- Explanation of Julia function(s) to be written (if applicable) 
    - Parameters/compatible inputs
    - Output returned by function
- Miscellaneous notes (if applicable)

#
## Q1: How different are the spike protein sequences from SARS-CoV2 and SARS-CoV?
#
### **Analysis Ideation**
According to Jason McLellan, a researcher at the University of Texas at Austin who studied the SARS-CoV2 spike protein, “SARS-CoV-2 binds ACE2 [FN: ACE2 is the angiotensin converting enzyme located in the intestines, kidney, testis, gallbladder, and heart that accounts for SARS-CoV-2’s pathological mechanisms] more strongly than does the virus that caused the severe acute respiratory syndrome outbreak in 2003” [(Researchers in China report structure of the novel coronavirus bound to its human target)][1]. 

[1]: https://cen.acs.org/biological-chemistry/biochemistry/Researchers-in-China-report-structure-of-the-novel-coronavirus-bound-to-its-human-target/98/web/2020/03
Because a protein’s amino acid sequencing determines its tertiary/quaternary structures which consequently determines its interactions with other proteins, this implies that the SARS-CoV2 and SARS-CoV spike proteins differ considerably in their sequencing. However, because these two viruses target the same host protein, ACE2, we also know that their sequences will not drastically different. I am curious to find out the degree of similarity between the spike protein sequences of these two viruses.

### **Sub-questions of Analysis**
**1. How well do these two protein sequences align?**
- **Julia function(s) utilized:** 
    - Needleman-Wunsch or Smith-Waterman alignment algorithm
        - Inputs: SARS-CoV2 spike protein sequence, SARS-CoV spike protein sequence
- **Method(s) for data visualization:**
    - a simplified version of a **[bioinformatics dot-plot][2]**, where the SARS-CoV2 spike protein sequence lies on the x-axis and the SARS-CoV spike protein sequence lies on the y-axis with the extreme N-terminal residues near the origin, and a dot is plotted for every position with matching residues; as a result, consecutively matching residues (i.e. matching subsequences) will be displayed as diagonal lines pointing from bottom left corner to upper right
    
    [2]: https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics)
**2. How many positions have (mis)matching residues?**
- **Julia function(s) utilized:**
    - `mis_matchData(alignedSeq)`
        - Input: output from NW or SW alignment algorithm
- **Method(s) for data visualization:**
    - A **bar chart** with two columns displaying the **total number** of matching (column 1) and mismatching (column 2) amino acid residues between SARS-CoV2 and SARS-CoV
    - A **pie chart** (with 2 slices) displaying the **relative percentages** of matching and mismatching amino acid residues between SARS-CoV2 and SARS-CoV

**3. What subsequences are (mis)matches and what are their lengths?**
- **Julia function(s) utilized:**
    - `mis_matchSubsequences(alignedSeq)`
        - Input: output from NW or SW alignment algorithm
- **Method(s) for data visualization:**
    - If the two spike protein sequences are extremely **similar** with **few mismatching sequences**, a **table** displaying the *positions* (e.g. 4-5), the *mismatching subsequences* and the *strain of virus they came from* (e.g. SARS-CoV2 - “LV”, SARS-CoV - “--”), and their *length* (e.g. 2 amino 
    acids)
    - If the two spike protein sequences are extremely **different** with **few matching subsequences**, then vice versa 

### **Julia functions to be written**
- `mis_matchData(alignedSeq)`
    - Compatible Input: takes a tuple of 2 aligned sequences as `Strings`
    - Compatible Output: a dictionary of `3 key-value` pairs, where the keys are the `Strings` `“Match#”`, `“Mismatch#”`, `“Mis/Match%”` and corresponding values are an Int of the number of total positions with **matching** residues, an Int of the number of total positions with **mismatching** residues, and a `Tuple` containing the **% of matches** and the **% of mismatches**, in this order. 
- `mis_matchSubsequences(alignedSeq)`
    - Compatible Input: takes a tuple of 2 aligned sequences, where each sequence is of type `String`
    - Compatible Output: a dictionary with an indeterminate number of key-value pairs where 
        - each value is an `Array` of **matching/mismatching protein subsequence(s)** as `Strings`
        - each key is a **3-element** `Tuple` with **1)** whether it is a “Match” or “Mismatch”, **2)** the corresponding position(s) of the subsequence as a string, and **3)** the length of this subsequence as an `Int`
    
    > Note 1: If the residue subsequences are a match, the value Array will contain 1 String, if they are a mismatch, the Array will contain 2 Strings

    > Note 2: for example, given an `alignedSeq` of (MES`LV`PGFNEKTH`V`QLS, MES`--`PGFNEKTH`A`QLS), `mis_matchSubsequences(alignedSeq)` returns {(“Match”, “1-3”, 3)=>[“MES”]; (“Mismatch”, “4-5”, 2)=>[“LV”, “--”]; (“Match”, “6-13”, 8)=>[“PGFNEKTH”]; (“Mismatch”, “14”, 1)=>[“V”, “A”]; (“Match”, “15-17”, 3)=>[“QLS”]}

#
## Q2: Is there a correlation between the distance between the location of where a SARS-CoV2 sample was collected and the location of the earliest reported COVID case in the U.S. (Seattle, WA) and the number of genomic residue mismatches? If so, is the correlation positive or negative? (i.e. positive = the farther away from Seattle, WA a sample is collected, the more different the SARS-CoV2 genomes are)
#
> Note 1: NCBI does not label samples with the exact location (e.g. City, State) of where they were collected, so the distance portion of this analysis will be the distance in miles between Seattle, WA, the location of the earliest reported COVID case in the U.S., and State A as displayed on Google maps

> Note 2: A single sample of SARS-CoV2 genome will be selected from every State, and each individual State sample will be compared with that of the earliest reported U.S. case. To minimize confounding factors, the individual State samples will be selected on the same day, which will be determined by a random number generator. Which sample is selected from each State’s database will also be determined by a random number generator.

- **Julia function(s) utilized:**
    - `Unique Kmers function`
        > Will need to run this function once for every State sample and once for the Seattle, WA sample 
        - Input (x 50): SARS-CoV2 genomic sequence from State A, B, C, D, E, ...
        - Input: the first SARS-CoV2 genomic sequence ever collected in Seattle, WA on January 21, 2020 
    - `Kmer set distance function`
        > will need to run this function 50 times to generate 50 distance metrics between each State sample and the Seattle, WA sample
        - Inputs: kmer sets from the SARS-CoV2 genome collected from a single sample in State A and the first sample ever collected in Seattle, WA on January 21, 2020 
- **Method(s) for data visualization:**
    - A **statistical dot plot** of **50 data pairs**, where each data pair consists of 
        1. the distance between the State of sample collection and Seattle, WA plotted on the Y-axis
        2. the distance metric of these SARS-CoV2 genomes plotted on the X-axis
    - A **linear regression** will be run with these **50 data pairs**, and the **line of best fit** resulting from this linear regression will be plotted on the statistical dot plot mentioned above




