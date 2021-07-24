# Final Project: Analysis 1

## Analysis Question: How different are the SARS-CoV and SARS-CoV2 spike / surface glycoprotein sequences?
* * *
## Goals & Hypothesis:
The goals of the following analysis is to answer the following questions:
1. How well do these two  sequences align?
2. How many positions have (mis)matching residues? 
3. In the aligned sequences, which subsequences are matches and at which positions do they occur?

Considering that SARS-CoV and SARS-CoV2 belong to the same genus and bind to the same receptor, ACE2, in hosts,
I hypothesize that their spike / surface glycoprotein sequences do not differ substantially. Because they have similar chemical properties and mechanisms of action,
most of their spike / surface glycoproteins should be conserved, especially regions that make up active sites.
* * *
## Methods Overview:
1. Parse the data
2. Align the sequences using the Needleman-Wunsch algorithm
3. Visualize this alignment with a bioinformatics dot-plot
4. Determine the counts and relative proportions of positions with matching and mismatching residues in the aligned sequences and visualize this with a bar chart and pie graph, respectively
5. Determine any matching subsequences longer than 20 residues in the aligned glycoprotein sequences
* * *
### 1. Setting Up
```julia
using BioinformaticsBISC195

CoV_spike_sequence = parse_fasta("data/SARS-cov-spike-glycoprotein.fasta"; DNA=false)
CoV2_surface_sequence = parse_fasta("data/SARS-cov2-surface-glycoprotein.fasta"; DNA=false)
```
**CoV_spike_sequence** =

```
  (Any["YP_009825051.1 |Severe acute respiratory syndrome-related coronavirus|spike glycoprotein|Canada: Toronto|"], Any  ["MFIFLLFLTLTSGSDLDRCTTFDDVQAPNYTQHTSSMRGVYYPDEIFRSDTLYLTQDLFLPFYSNVTGFHTINHTFGNPVIPFKDGIYFAATEKSNVVRGWVFGSTMNNKSQSVIIINNSTNVVIRACNFELCDNPFFAVSKPMGTQTHTMIFDNAFNCTFEYISDAFSLDVSEKSGNFKHLREF  VFKNKDGFLYVYKGYQPIDVVRDLPSGFNTLKPIFKLPLGINITNFRAILTAFSPAQDIWGTSAAAYFVGYLKPTTFMLKYDENGTITDAVDCSQNPLAELKCSVKSFEIDKGIYQTSNFRVVPSGDVVRFPNITNLCPFGEVFNATKFPSVYAWERKKISNCVADYSVLYNSTFFSTFKCYGVSAT  KLNDLCFSNVYADSFVVKGDDVRQIAPGQTGVIADYNYKLPDDFMGCVLAWNTRNIDATSTGNYNYKYRYLRHGKLRPFERDISNVPFSPDGKPCTPPALNCYWPLNDYGFYTTTGIGYQPYRVVVLSFELLNAPATVCGPKLSTDLIKNQCVNFNFNGLTGTGVLTPSSKRFQPFQQFGRDVSDFT  DSVRDPKTSEILDISPCAFGGVSVITPGTNASSEVAVLYQDVNCTDVSTAIHADQLTPAWRIYSTGNNVFQTQAGCLIGAEHVDTSYECDIPIGAGICASYHTVSLLRSTSQKSIVAYTMSLGADSSIAYSNNTIAIPTNFSISITTEVMPVSMAKTSVDCNMYICGDSTECANLLLQYGSFCTQLN  RALSGIAAEQDRNTREVFAQVKQMYKTPTLKYFGGFNFSQILPDPLKPTKRSFIEDLLFNKVTLADAGFMKQYGECLGDINARDLICAQKFNGLTVLPPLLTDDMIAAYTAALVSGTATAGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKQIANQFNKAISQIQESLTTTSTALGKLQDV  VNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQAAPHGVVFLHVTYVPSQERNFTTAPAICHEGKAYFPREGVFVFNGTSWFITQRNFFSPQIITTDNTFVSGNCDVVIGIINNTVY  DPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYVWLGFIAGLIAIVMVTILLCCMTSCCSCLKGACSCGSCCKFDEDDSEPVLKGVKLHYT"])
  ```
**CoV2_surface_sequence** =
  ```
  (Any["YP_009724390.1 |Severe acute respiratory syndrome-related coronavirus|surface glycoprotein|China|"], Any["MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"])
  ```
  ***
### 2. Aligning the Two Glycoprotein Sequences
```julia
  aligned = nwalign(CoV_spike_sequence[2][1], CoV2_surface_sequence[2][1])
```
**aligned** =
  ```
  ("MFIFL-LFLTLTSGSDLDRCTTFDDVQAPN-YTQHTSSMRGVYYPDEIFRSDTLYLTQDLFLPFYSNVTGFHTI-----NHT--FGNPVIPFKDGIYFAATEKSNVVRGWVFGSTMNNKSQSVIIINNSTNVVIRACNFELC-DNPFFAVS--KPMGTQTHTMIFDN---AFNCTFEYISDAFSL-DVSE-KSGNFKHLREFVFKNKDG-FL-YVYKGYQPIDVVRDLPSGFNTLKPIFKLPLGINITNFRAILTAF--S---PAQDIW-G-TS-AAAYFVGYLKPTTFMLKYDENGTITDAVDCSQNPLAELKCSVKSF--EIDKGIYQTSNFRVVP--SGDVVRFPNITNLCPFGEVFNATKFPSVYAWERKKISNCVADYSVLYNSTFFSTFKCYGVSATKLNDLCFSNVYADSFVVKGDDVRQIAPGQTGVIADYNYKLPDDFMGCVLAWNTRNIDATST-G-NYNYKYRYL-RHGK--LRPFERDISNVPFSPDGK-PCTPPAL-NCYWPLNDYGFYTTTGIGYQPYRVVVLSFELLNAPATVCGPKLSTDLIKNQCVNFNFNGLTGTGVLTPSSKRFQPFQQFGRDVSDFTDSVRDPKTSEILDISPCAFGGVSVITPGTNASSEVAVLYQDVNCTDVSTAIHADQLTPAWRIYSTGNNVFQTQAGCLIGAEHVDTSYECDIPIGAGICASYHTV--SLLR--ST-SQKSIVAYTMSLGADSSIAYSNNTIAIPTNFSISITTEVMPVSMAKTSVDCNMYICGDSTECANLLLQYGSFCTQLNRALSGIAAEQDRNTREVFAQVKQMYKTPTLKYFGGFNFSQILPDPLKPTKRSFIEDLLFNKVTLADAGFMKQYGECLGDINARDLICAQKFNGLTVLPPLLTDDMIAAYTAALVSGTATAGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKQIANQFNKAISQIQESLTTT-STALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQAAPHGVVFLHVTYVPSQERNFTTAPAICHEGKAYFPREGVFVFNGTSWFITQRNFFSPQIITTDNTFVSGNCDVVIGIINNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYVWLGFIAGLIAIVMVTILLCCMTSCCSCLKGACSCGSCCKFDEDDSEPVLKGVKLHYT", "MFVFLVL-LPLVS-SQCVNLTTRT--QLPPAYTN--SFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCND-PFLGVYYHKNNKSWMESE-FRVYSSANNCTFEYVSQPF-LMDL-EGKQGNFKNLREFVFKNIDGYFKIYS-KHT-PINLVRDLPQGFSALEPLVDLPIGINITRFQTLL-ALHRSYLTPG-DSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVE--KGIYQTSNFRVQPTESI-V-RFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLD--SKVGGNYNYLYR-LFR--KSNLKPFERDISTEIYQA-GSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQ-SIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTAS-ALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT")
  ```
***
### 3. Creating a Bioinformatics Dot-Plot

The CoV spike glycoprotein's residues line the x-axis, with residues at the extreme N-terminus towards the origin.

Each residue's one letter abbreviation is labeled with its position in front (is informative, but also because Plots.jl requires x (and y) points to be sorted in ascending order).
```julia
x = [string(i, aligned[1][i]) for i in 1:length(aligned[1])]
```
The CoV2 surface glycoprotein's residues line the y-axis, with residues at the extreme N-terminus towards the origin.

Each residue's one letter abbreviation is labeled with its position in front (is informative, but also because Plots.jl requires y (and x) points to be sorted in ascending order).
```julia
y = [string(i, aligned[2][i]) for i in 1:length(aligned[2])] 
```
The following function generates a matrix that serves as the "color-coding" for the heatmap.

The cell at matrix pos (1, 1) corresponds to the grid in the uppermost left corner of the heatmap.

A number of 0.0 in the matrix corresponds to a blank grid while a number of 1.0 corresponds to a grid filled in black.
```julia
z = alignmentHeatmap(aligned[1], aligned[2])
```
***
Using Plots.jl, we can generate a heatmap where
  - the origin is in the upper left corner (i.e. imagine a coordinate plane where only the fourth quadrant is visible)
  - two **aligned** sequences line the x and y axes
  - residues of the CoV spike glycoprotein are lined on the y-axis, with the residue at the extreme N-terminus at the origin
  - residues of the CoV2 surface glycoprotein are lined on the x-axis, with the residue at the extreme N-terminus at the origin
  - grids that correspond to identical residues are filled in _black_ while those that do not are left _blank_

This heatmap achieves the same goal in data visualization as a bioinformatics dot-plot, 
  which is used to compare two biologically sequences already aligned and identify regions of similarity.

If the two biological sequences are perfectly identical, 
  the diagonal line stretching the origin to the lower right corner of the plot should be clearly displayed and completely continous.

```julia
using Plots
heatmap(x, y, z, xmirror = true, yflip = true, color = cgrad([:white, :black]), xlabel = "Residues of CoV2 Surface Glycoprotein", ylabel = "Residues of CoV Spike Glycoprotein", title = "Bioinformatics Dot Plot\nComparing Different Coronaviruses' S Glycoprotein Sequence\n ", titlefontsize = 12, labelpadding = 2.0)
```

From the heatmap/bioinformatics dot-plot above,
  we can see that the CoV spike glycoprotein and CoV2 surface glycoprotein are extremely similar in sequencing,
  as the diagonal line going from the origin to the bottom right corner
  is easily discernible with only a few faint and broken regions in the middle.



***
### 4. Comparing the relative number of positions with matching and mistaching residues and visualizing the data

Given that the primary sequence of a
  polypeptide plays a crucial role in
  determining its
  three-dimensional structure and, therefore chemical properties (e.g. hydrophobicity, compatible substrates, binding affinities), determining the counts and relative proportions of positions with (mis)matching residues in the aligned glycoprotein sequences can not only help us better understand how (dis)similar they are genetically, but also functionally. 

```julia
# Generating a dictionary with the counts and relative proportions of positions with (mis)matching residues
mis_match_stats = mis_matchSeq(aligned)

# Creating the necessary components for organizing the data above in tabular form with DataFrames.jl
labels = ["Match", "Mismatch"]
counts = [mis_match_stats["Match#"], mis_match_stats["Mismatch#"]]
percentages = [(string(mis_match_stats["Mis/Match%"][1]*100)[1:4] * "%"), (string(mis_match_stats["Mis/Match%"][2]*100)[1:4] * "%")]
ratios = [mis_match_stats["Mis/Match%"][1], mis_match_stats["Mis/Match%"][2]]

using DataFrames

df_residue_alignment = DataFrame(Category=labels, Counts=counts, Ratios = ratios, Percentages=percentages)
```
***
To visualize the **relative proportions** of positions with (mis)matching residues in the aligned glycoprotein sequences, we can make a **pie graph**.
```julia
using VegaLite

df_residue_alignment |> 
@vlplot(:arc, theta={:Counts, stack=true}, color= {"Category:n"}, view={stroke=nothing}, title= {text = "Relative Proportions of Matching and Matching Residues in Protein Sequences", subtitle = "CoV Spike Protein and CoV2 Surface Protein"}) + 
@vlplot(mark={:arc, outerRadius=80}) + 
@vlplot(mark={:text, radius=95}, text="Percentages:n")
```
***
To visualize the **counts** of positions with (mis)matching residues in the aligned glycoprotein sequences, we can make a **bar graph**.
```julia
using Plots 

bar([(labels[1], counts[1]), (labels[2], counts[2])], legend = false, bar_width = .5, xlabel = "Residue Alignment", ylabel = "# of Residues", title = "Counts of Matching and Mismatching Residue Alignment Of\nCoV Spike Glycoprotein and CoV2 Surface Glycoprotein\n ", titlefontsize = 12)
```
Looking at the data table, bar chart, and pie graph from above, we see that nearly 80% of the two glycoproteins' residues (nearly 1,000 out of the total ~1,300 residues) are identical at the same positions. This shows that the SARS-CoV and SARS-CoV2 spike / surface glycoproteins' primary sequences are predominantly the same, from which we can infer that the two glycoproteins are also similar in three-dimensional structure and therefore in functional properties. 
***
### 5. Determining Any Matching Subsequences Longer Than 20 Residues In the Aligned Glycoprotein Sequences

In Step #4, we looked at the two
  glycoprotein sequences' residues as individual units and discovered that nearly 80% of them are identical at the same positions. 

But, how many of these matching residues occur at consecutive positions? From the data generated in Step #4, we could conclude that the two sequences are highly conserved. However, this conclusion would be premature if we did not look at the number of matching residues that occur at consecutive positions, or, in other words, the lengths of the matching *sub*sequences.

If two sequences are highly conserved, they should have a substantial number of *long* (arbitrarily defined as 20 residues here) matching subsequences.

```julia
using BioinformaticsBISC195
subsequences_stats = mis_matchSubsequences(aligned)

count(item->(item.type=="match" && length(item.sequences[1])>20), subsequences_stats)

```
After running the script above, we can see that there are **6** *long* matching subsequences in the aligned glycoprotein sequences.

However, it is of even greater interest to determine the identities and positions of these subsequences, as these bits of information can 
- give us a preliminary understanding of the locations of highly-conserved regions of the SARS-CoV and SARS-CoV2 spike / surface glycoproteins 
- help us analyze the intermolecular interactions that may arise between the R-groups of the residues in these highly conserved regions, which may help us further analyze the overall chemical properties and functions of these regions holistically after the polypeptides are folded into their three-dimensional structures

```julia
# Creating the necessary components for organizing this data in tabular form with DataFrames.jl
begin 
    match_20_up_arr = []
    for subseq in subsequences_stats
        if subseq.type=="match" && length(subseq.sequences[1])>20
            push!(match_20_up_arr, (subseq.positions, subseq.sequences[1]))
        end
    end
end


using DataFrames

df_subsequence_match_20up = DataFrame(Positions = [i[1] for i in match_20_up_arr], Sequences = [i[2] for i in match_20_up_arr])
```
| Positions      | Sequences |
| ----------- | ----------- |
| 871:892      | ARDLICAQKFNGLTVLPPLLTD       |
| 910:946   | GWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQK        |
| 970:1080   | ALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRL
| |QSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQ        |
| 1138:1158   | PQIITTDNTFVSGNCDVVIGI        |
| 1160:1241   | NNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWY        |
| 1274:1299   | CSCGSCCKFDEDDSEPVLKGVKLHYT        |

From the table above, we can see that the majority of the highly-conserved subsequences occur towards the C-terminus in the aligned glycoprotein sequences (1299 residues long).

Further analysis is out of the scope of my knowledge, but I suspect that 
- these subsequences may be crucial to the chemical properties of the coronaviruses such that mutations here significantly impact the viruses' ability to enter host cells, survive, and reproduce
- these subsequences may constitute the active sites of the folded polypeptides
- the intermolecular interactions of the residues in these subsequences may be crucial in constructing or maintaining the protein's three-dimensional structure (e.g. if some regions mainly comprise of acidic amino acids and others basic amino acids such that the dipole-dipole attraction gives rise to the protein's three-dimensional structure)