<!-- TODO: Better name for this file, use headers / sections -->
```julia
using BioinformaticsBISC195

CoV_spike_sequence = parse_fasta("data/SARS-cov-spike-glycoprotein.fasta"; DNA=false)
CoV2_surface_sequence = parse_fasta("data/SARS-cov2-surface-glycoprotein.fasta"; DNA=false)

aligned = nwalign(CoV_spike_sequence[2][1], CoV2_surface_sequence[2][1])
```
- **CoV_spike_sequence** =

  ```
  (Any["YP_009825051.1 |Severe acute respiratory syndrome-related coronavirus|spike glycoprotein|Canada: Toronto|"], Any  ["MFIFLLFLTLTSGSDLDRCTTFDDVQAPNYTQHTSSMRGVYYPDEIFRSDTLYLTQDLFLPFYSNVTGFHTINHTFGNPVIPFKDGIYFAATEKSNVVRGWVFGSTMNNKSQSVIIINNSTNVVIRACNFELCDNPFFAVSKPMGTQTHTMIFDNAFNCTFEYISDAFSLDVSEKSGNFKHLREF  VFKNKDGFLYVYKGYQPIDVVRDLPSGFNTLKPIFKLPLGINITNFRAILTAFSPAQDIWGTSAAAYFVGYLKPTTFMLKYDENGTITDAVDCSQNPLAELKCSVKSFEIDKGIYQTSNFRVVPSGDVVRFPNITNLCPFGEVFNATKFPSVYAWERKKISNCVADYSVLYNSTFFSTFKCYGVSAT  KLNDLCFSNVYADSFVVKGDDVRQIAPGQTGVIADYNYKLPDDFMGCVLAWNTRNIDATSTGNYNYKYRYLRHGKLRPFERDISNVPFSPDGKPCTPPALNCYWPLNDYGFYTTTGIGYQPYRVVVLSFELLNAPATVCGPKLSTDLIKNQCVNFNFNGLTGTGVLTPSSKRFQPFQQFGRDVSDFT  DSVRDPKTSEILDISPCAFGGVSVITPGTNASSEVAVLYQDVNCTDVSTAIHADQLTPAWRIYSTGNNVFQTQAGCLIGAEHVDTSYECDIPIGAGICASYHTVSLLRSTSQKSIVAYTMSLGADSSIAYSNNTIAIPTNFSISITTEVMPVSMAKTSVDCNMYICGDSTECANLLLQYGSFCTQLN  RALSGIAAEQDRNTREVFAQVKQMYKTPTLKYFGGFNFSQILPDPLKPTKRSFIEDLLFNKVTLADAGFMKQYGECLGDINARDLICAQKFNGLTVLPPLLTDDMIAAYTAALVSGTATAGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKQIANQFNKAISQIQESLTTTSTALGKLQDV  VNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQAAPHGVVFLHVTYVPSQERNFTTAPAICHEGKAYFPREGVFVFNGTSWFITQRNFFSPQIITTDNTFVSGNCDVVIGIINNTVY  DPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYVWLGFIAGLIAIVMVTILLCCMTSCCSCLKGACSCGSCCKFDEDDSEPVLKGVKLHYT"])
  ```
- **CoV2_surface_sequence** =
  ```
  (Any["YP_009724390.1 |Severe acute respiratory syndrome-related coronavirus|surface glycoprotein|China|"], Any["MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"])
  ```
- **aligned** =
  ```
  ("MFIFL-LFLTLTSGSDLDRCTTFDDVQAPN-YTQHTSSMRGVYYPDEIFRSDTLYLTQDLFLPFYSNVTGFHTI-----NHT--FGNPVIPFKDGIYFAATEKSNVVRGWVFGSTMNNKSQSVIIINNSTNVVIRACNFELC-DNPFFAVS--KPMGTQTHTMIFDN---AFNCTFEYISDAFSL-DVSE-KSGNFKHLREFVFKNKDG-FL-YVYKGYQPIDVVRDLPSGFNTLKPIFKLPLGINITNFRAILTAF--S---PAQDIW-G-TS-AAAYFVGYLKPTTFMLKYDENGTITDAVDCSQNPLAELKCSVKSF--EIDKGIYQTSNFRVVP--SGDVVRFPNITNLCPFGEVFNATKFPSVYAWERKKISNCVADYSVLYNSTFFSTFKCYGVSATKLNDLCFSNVYADSFVVKGDDVRQIAPGQTGVIADYNYKLPDDFMGCVLAWNTRNIDATST-G-NYNYKYRYL-RHGK--LRPFERDISNVPFSPDGK-PCTPPAL-NCYWPLNDYGFYTTTGIGYQPYRVVVLSFELLNAPATVCGPKLSTDLIKNQCVNFNFNGLTGTGVLTPSSKRFQPFQQFGRDVSDFTDSVRDPKTSEILDISPCAFGGVSVITPGTNASSEVAVLYQDVNCTDVSTAIHADQLTPAWRIYSTGNNVFQTQAGCLIGAEHVDTSYECDIPIGAGICASYHTV--SLLR--ST-SQKSIVAYTMSLGADSSIAYSNNTIAIPTNFSISITTEVMPVSMAKTSVDCNMYICGDSTECANLLLQYGSFCTQLNRALSGIAAEQDRNTREVFAQVKQMYKTPTLKYFGGFNFSQILPDPLKPTKRSFIEDLLFNKVTLADAGFMKQYGECLGDINARDLICAQKFNGLTVLPPLLTDDMIAAYTAALVSGTATAGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKQIANQFNKAISQIQESLTTT-STALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQAAPHGVVFLHVTYVPSQERNFTTAPAICHEGKAYFPREGVFVFNGTSWFITQRNFFSPQIITTDNTFVSGNCDVVIGIINNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYVWLGFIAGLIAIVMVTILLCCMTSCCSCLKGACSCGSCCKFDEDDSEPVLKGVKLHYT", "MFVFLVL-LPLVS-SQCVNLTTRT--QLPPAYTN--SFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCND-PFLGVYYHKNNKSWMESE-FRVYSSANNCTFEYVSQPF-LMDL-EGKQGNFKNLREFVFKNIDGYFKIYS-KHT-PINLVRDLPQGFSALEPLVDLPIGINITRFQTLL-ALHRSYLTPG-DSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVE--KGIYQTSNFRVQPTESI-V-RFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLD--SKVGGNYNYLYR-LFR--KSNLKPFERDISTEIYQA-GSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQ-SIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTAS-ALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT")
  ```


- CoV spike glycoprotein's residues line the x-axis, with residues at the extreme N-terminus towards the origin
- Each residue's one letter abbreviation is labeled with its position in front (is informative, but also because Plots.jl requires x (and y) points to be sorted in ascending order)
```julia
x = [string(i, aligned[1][i]) for i in 1:length(aligned[1])]
```
- CoV2 surface glycoprotein's residues line the y-axis, with residues at the extreme N-terminus towards the origin
- Each residue's one letter abbreviation is labeled with its position in front (is informative, but also because Plots.jl requires y (and x) points to be sorted in ascending order)
```julia
y = [string(i, aligned[2][i]) for i in 1:length(aligned[2])] 
```
- The following function generates a matrix that serves as the "color-coding" for the heatmap
- The cell at matrix pos (1, 1) corresponds to the grid in the uppermost left corner of the heatmap
- A number of 0.0 in the matrix corresponds to a blank grid while a number of 1.0 corresponds to a grid filled in black
```julia
z = alignmentHeatmap(aligned[1], aligned[2])
```
- Using Plots.jl, we can generate a heatmap where
    - the origin is in the upper left corner (i.e. imagine a coordinate plane where only the fourth quadrant is visible)
    - two **aligned** sequences line the x and y axes
        - residues of the CoV spike glycoprotein are lined on the y-axis, with the residue at the extreme N-terminus at the origin
        - residues of the CoV2 surface glycoprotein are lined on the x-axis, with the residue at the extreme N-terminus at the origin
    - grids that correspond to identical residues are filled in _black_ while those that do not are left _blank_
- This heatmap achieves the same goal in data visualization as a bioinformatics dot-plot, which is used to compare two biologically sequences already aligned and identify regions of similarity
    - If the two biological sequences are perfectly identical, the diagonal line stretching the origin to the lower right corner of the plot should be clearly displayed and completely continous
```julia
using Plots


heatmap(x, y, z, xmirror = true, yflip = true, color = cgrad([:white, :black]),
    xlabel = "Residues of CoV2 Surface Glycoprotein",
    ylabel = "Residues of CoV Spike Glycoprotein")
```
- From the heatmap/bioinformatics dot-plot above,
  we can see that the CoV spike glycoprotein and CoV2 surface glycoprotein are extremely similar in sequencing,
  as the diagonal line going from the origin to the bottom right corner
  is easily discernible with only a few faint and broken regions in the middle

heatmap(x, y, z, xmirror = true, yflip = true, color = cgrad([:white, :black]), xlabel = "Residues of CoV2 Surface Glycoprotein", ylabel = "Residues of CoV Spike Glycoprotein", title = "Bioinformatics Dot Plot\nComparing Different Coronaviruses' S Glycoprotein Sequence\n ", titlefontsize = 12, labelpadding = 2.0)
```
- From the heatmap/bioinformatics dot-plot above, we can see that the CoV spike glycoprotein and CoV2 surface glycoprotein are extremely similar in sequencing, as the diagonal line going from the origin to the bottom right corner is easily discernible with only a few faint and broken regions in the middle

#
```julia
mis_match_stats = mis_matchSeq(aligned)
Dict{String, Any} with 3 entries:
  "Mis/Match%" => (0.767513, 0.232487)
  "Mismatch#"  => 302
  "Match#"     => 997
```
```julia
labels = ["Match", "Mismatch"]
counts = [mis_match_stats["Match#"], mis_match_stats["Mismatch#"]]
percentages = [(string(mis_match_stats["Mis/Match%"][1]*100)[1:4] * "%"), (string(mis_match_stats["Mis/Match%"][2]*100)[1:4] * "%")]
ratios = [mis_match_stats["Mis/Match%"][1], mis_match_stats["Mis/Match%"][2]]
```
```julia
using DataFrames

df_residue_alignment = DataFrame(Category=labels, Counts=counts, Ratios = ratios, Percentages=percentages)
2×4 DataFrame
 Row │ Category  Counts  Ratios    Percentages 
     │ String    Int64   Float64   String      
─────┼─────────────────────────────────────────
   1 │ Match        997  0.767513  76.7%
   2 │ Mismatch     302  0.232487  23.2%
```
```julia
df_residue_alignment |> 
@vlplot(:arc, theta={:Counts, stack=true}, color= {"Category:n"}, view={stroke=nothing}, title= {text = "Relative Proportions of Matching and Matching Residues in Protein Sequences", subtitle = "CoV Spike Protein and CoV2 Surface Protein"}) + 
@vlplot(mark={:arc, outerRadius=80}) + 
@vlplot(mark={:text, radius=95}, text="Percentages:n")
```
```julia
bar([(labels[1], counts[1]), (labels[2], counts[2])], legend = false, bar_width = .5, xlabel = "Residue Alignment", ylabel = "# of Residues", title = "Counts of Matching and Mismatching Residue Alignment Of\nCoV Spike Glycoprotein and CoV2 Surface Glycoprotein\n ", titlefontsize = 12)
```

```julia
subsequences_stats = mis_matchSubsequences(aligned)
395-element Vector{Any}:
 AlignmentSubsequences("match", 1:2, ["MF"])
 AlignmentSubsequences("mismatch", 3:3, ["I", "V"])
 AlignmentSubsequences("match", 4:5, ["FL"])
 AlignmentSubsequences("mismatch", 6:6, ["-", "V"])
 AlignmentSubsequences("match", 7:7, ["L"])
 AlignmentSubsequences("mismatch", 8:8, ["F", "-"])
 ⋮
 AlignmentSubsequences("match", 1243:1258, ["WLGFIAGLIAIVMVTI"])
 AlignmentSubsequences("mismatch", 1259:1259, ["L", "M"])
 AlignmentSubsequences("match", 1260:1272, ["LCCMTSCCSCLKG"])
 AlignmentSubsequences("mismatch", 1273:1273, ["A", "C"])
 AlignmentSubsequences("match", 1274:1299, ["CSCGSCCKFDEDDSEPVLKGVKLHYT"])
```

```julia
count(item->(item.type=="match" && length(item.sequences[1])>20), subsequences_stats)
6
```


```julia
begin 
    match_20_up_arr = []
    for subseq in subsequences_stats
        if subseq.type=="match" && length(subseq.sequences[1])>20
            push!(match_20_up_arr, (subseq.positions, subseq.sequences[1]))
            @info subseq.positions
            @info subseq.sequences[1]
        end
    end
    return match_20_up_arr
end

6-element Vector{Any}:
 (871:892, "ARDLICAQKFNGLTVLPPLLTD")
 (910:946, "GWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQK")
 (970:1080, "ALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQ")
 (1138:1158, "PQIITTDNTFVSGNCDVVIGI")
 (1160:1241, "NNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWY")
 (1274:1299, "CSCGSCCKFDEDDSEPVLKGVKLHYT")
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
```julia
df_subsequence_match_20up = DataFrame(Positions = [i[1] for i in match_20_up_arr], Sequences = [i[2] for i in match_20_up_arr])
6×2 DataFrame
 Row │ Positions  Sequences                         
     │ UnitRang…  String                            
─────┼──────────────────────────────────────────────
   1 │ 871:892    ARDLICAQKFNGLTVLPPLLTD
   2 │ 910:946    GWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQK
   3 │ 970:1080  ALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQ  
   4 │ 1138:1158  PQIITTDNTFVSGNCDVVIGI
   5 │ 1160:1241  NNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWY
   6 │ 1274:1299  CSCGSCCKFDEDDSEPVLKGVKLHYT
```
