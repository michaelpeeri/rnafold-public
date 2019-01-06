
## Supplementary Figure S1. All species ΔLFE on taxonomic tree

![](alltaxa.300nt.slice_0_0.png)
![](alltaxa.300nt.slice_1_0.png)
![](alltaxa.300nt.slice_2_0.png)
![](alltaxa.300nt.slice_3_0.png)
![](alltaxa.300nt.slice_4_0.png)
![](alltaxa.300nt.slice_5_0.png)
![](alltaxa.300nt.slice_6_0.png)
![](alltaxa.300nt.slice_7_0.png)
![](alltaxa.300nt.slice_8_0.png)
![](alltaxa.300nt.slice_9_0.png)

ΔLFE profiles calculated using the CDS-wide randomization
for individual species arranged by NCBI taxonomy. The
ΔLFE profiles shown are for positions 0-300nt relative to
CDS start (left) and CDS end (right). The numbers of species
included in each group is shown to the left of the group name.


## Supplementary Figure S2. Comparison between ΔLFE calculated using CDS-wide and position-specific ("vertical") randomization

![](correlation_between_series.vertical.combined.png){ width=60% }
![](heatmap_profile_legend.png){ width=20% } ΔLFE (kcal/mol/window)

Position-specific randomization (maintaining the encoded AA sequences as well as the codon frequency in each position (across all CDSs belonging to the same species) yields qualitatively similary results to the CDS-wide randomization used throughout the rest of this paper. This supports the conclusion that the observed ΔLFE profiles are not merely a result of position-dependent biases in codon composition.
**A** Correlation between ΔLFE calculated using "CDS-wide" and "position-specific" randomizations (see methods), at each position relative to CDS start. Correlations were calculated for a random sample (*N*=23) of species.
**B** Comparison of individual mean ΔLFE profiles calculated using "CDS-wide" (LFE-0) and "position-specific" (LFE-1) randomizations.

## Supplementary Figure S3. Traits correlogram (correlation vs. phylogenetic distance)

![](Correlogram.1.pdf.png)

Correlation (expressed using Moran's I coefficient) between the values of different traits, for pairs of species of different phylogenetic distances. Genomic-GC% is positively correlated at short distances. ΔLFE values (at different positions relative to CDS start) are more strongly correlated than genomic-GC% at most phylogenetic distances, but less correlated than genome sizes. Confidence intervals represent 95% confidence calculated using 500 bootstrap samples. The 'Random' trait is a normally-distributed uncorrelated variable.

## Supplementary Figure S4. Local CUB vs. Local ΔLFE

![](cub_corrs_all_hist_only.pdf.png){ width=30% }

Spearman correlations between the ΔLFE profile and the corresponding CUB profiles show no direct correspondence, indicating the ΔLFE profiles are not simply a side-effect of direct selection operating on CUB in different CDS regions. CUB measures were calculated for the sequences contained in the same 40nt windows, starting at positions 0-300nt relative to CDS start, with all the sequences for each species concatenated, for a random sample of *N*=256 species.

From top to bottom, Nc (Effective Number of Codons), CAI (Codon Adaptation Index), Fop (Frequency of Optimal Codons), GC% (GC-content).


## Supplementary Figure S5. Unsupervised discovery of ΔLFE profile regions

**A**
![](pca_profiles_combined_begin.png){ width=30% }
**B**
![](pca_profiles_combined_end.png){ width=30% }  
ΔLFE
![](heatmap_profile_legend.png){ width=20% }  
**C**
![](pca_profiles_combined_angle_1_2.png){ width=30% }  

Principal Component Analysis (PCA) of the ΔLFE profiles uncovers two components, with different relative weights for the CDS-edge and mid-CDS regions. **A** PCA plot for ΔLFE profiles at positions 0-300nt relative to CDS start (represented as vectors of length 31), shown by plotting each ΔLFE profile in its position in PCA space (with 2 dimensions), with overlapping profiles hidden to avoid clutter. The density of profiles in each region is illustrated using shading and the marginal distributions are shown on the axes. Loading vectors for positions 0nt and 250nt (relative to CDS start) are shown. To verify this analysis is robust, bootstrapping using 1000 repeats was used to measure the following values:
RSD1 - Relative standard-deviation (SD/mean) for the angle between the loading vectors shown (i.e., those for ΔLFE profile positions 0nt and 250nt). Distribution of angles shown in **C**.
RSD2 - Relative standard-deviation (SD/mean) for the explained variance of PC1.
**B** PCA plot for ΔLFE profiles at positions 0-300nt relative to CDS end (created using the same method as **A**).
**C**
Distribution of angles between shown loading vectors (i.e., those for ΔLFE profile positions 0nt and 250nt) using 1000 bootstrap samples. The distribution mean is 2.08 radians ($119^{\circ}$) and the relative standard deviation (also shown as RSD1 on **A**) is 1.4%.

This procedure was repeated for all species and for each domain individially (see also Fig. 4D). In each case, the first two PCs explain >80% of the variation. *TODO* Give max values for RSD1 and RSD2.
The loading vectors for positions 0nt and 250nt are not parralel nor orthogonal (and this is robust to sampling and persists in smaller groups, see Figure 4D), indicating an some level of dependence between the two positions (also indicated in Fig. 3E, Supplementary Fig. S6).

## Supplementary Figure S6. Autocorrelation between ΔLFE profile regions

**A**
![](Correlations_begin_300nt.png){ width=50% }
**B**
![](Correlations_table_begin_300nt.png){ width=35% }

**A**
Autocorrelation for ΔLFE between positions relative to CDS start.
Above main diagonal - Pearson's correation. Below main diagonal - coefficient of determination ($R^2$) for GLS regression.
Values for positions a-h indicated on **B**.
Significant positions (*p*-value<0.01) indicated by white dots.
**B**
Numerical values (a-d - $R^2$, e-h - Pearson's-*r*) and *p*-values for positions marked in **A**.  
This supports the robustness of the values in Fig. 3E.

## Supplementary Figure S7. Genomic-GC% correlation profile does not come from components

![](R_comparison_dLFE_native_shuffled-nup.pdf.png)

Coefficient of determination ($R^2$) with regression direction for GLS regression of the specified trait with ΔLFE and its components (ΔLFE - red; native LFE - green; randomized LFE - blue), at different positions relative to CDS start. The observed correlation between each trait and ΔLFE is not observed with the individial components.

## Supplementary Figure S8. Trait correlations in taxonomic subgroups

**A**
![](regressionByTaxgroup_GenomicGC_range_1_2.pdf.png){ width=60% }  
**B**
![](regressionByTaxgroup_GenomicENc.prime_range_1_2_reduced.svg.png){ width=30% }
**C**
TODO
**D**  
TODO

Coefficient of determination ($R^2$) and regression direction for GLS regression between genomic-GC% and mean ΔLFE in different taxonomic subgroups, for two regions relative to CDS-start. Top bar - 0-20nt; Bottom bar - 70-300nt. Sign of regression slope is indicated by color - Red - positive (reinforcing) effect; Blue - negative (compensating) effect. Significant results (*p*-value<0.01) are indicated by color intensity and marked with a '$\ast$'. Included taxonomic groups have 9 or more species in the dataset.

**A** Genomic GC
**B** Genomic ENc'
**C** Optimum Temperature
**D** TODO

## Supplementary Figure S9. Endosymbionts have weaker ΔLFE in the mid-CDS region

**A**
![](Endosymbiont_comparison.pdf.png){ width=40% }
**B**
![](Endosymbiont_comparison_end.pdf.png){ width=40% }

Comparison between mean abs(ΔLFE) for the region 150-300nt relative to CDS start, for known endosymbionts (*N*=11 included in the tree) (see Supplementary Table S1) and species not known as endosymbionts (*N*=325). ΔLFE values are normalized to have mean=0. Known endosymbionts have significantly higher mean ΔLFE in this region (GLS, *p*-value=9e-4, $R^2$=0.04). **B** Comparison between mean abs(ΔLFE) for the region 150-300nt relative to CDS end, for known endosymbionts (using the same method as **A**).


## Supplementary Figure S10. Genomic-ENc' correlates with ΔLFE magnitude, not shape

**A**
![](R2comparison_GenomicENc.prime_profile_1_all.pdf.png){ width=40% }
**B**
![](R2comparison_GenomicENc.prime_profile_2_all.pdf.png){ width=40% }  
**C**
![](R2comparison_GenomicGC_profile_1_all.modified.svg.png){ width=40% }
**D**
![](R2comparison_GenomicGC_profile_2_all.modified.svg.png){ width=40% }  

The observed correlation of ΔLFE with Genomic-ENc' (Figure 6) is due to correlation with the magnitude of the ΔLFE profile. When all profiles are normalized to have the same scale (by dividing the values of each profile by their standard deviation so the resulting profiles all have standard deviation 1), most of the correlation is removed (**A**,**B**). For comparison, the same procedure is followed for genomic-GC (**C**,**D**). Values represent coefficient of determination ($R^2$) for GLS regression of each trait (genomic-ENc' or genomic-GC%) vs. the normalized ΔLFE profile at different position relative to CDS edges, with the sign representing the regression coefficient. Regressions for different taxons are shown using different line colors and widths (black is for all species), and white dots show areas in which the regression is significant (*p*-value<0.01). The dashed red line represent $R^2$ for regression against the standard deviation for each ΔLFE profile (i.e., the scaling factor).
**A** Genomic-ENc' vs. ΔLFE, CDS start.
**B** Genomic-ENc' vs. ΔLFE, CDS end.
**C** Genomic-GC vs. ΔLFE, CDS start.
**D** Genomic-GC vs. ΔLFE, CDS end.


## Supplementary Figure S11. Range robustness for GLS regressions between ΔLFE and related traits

**A**
![](tree_phenotypes_regression.out.dLFE.by.taxgroups.78-nup.11.png){ width=30% }
![](tree_phenotypes_regression.out.dLFE.by.taxgroups.78-nup.12.png){ width=30% }
![](tree_phenotypes_regression.out.dLFE.by.taxgroups.78-nup.13.png){ width=30% }
![](tree_phenotypes_regression.out.dLFE.by.taxgroups.78-nup.14.png){ width=30% }
![](tree_phenotypes_regression.out.dLFE.by.taxgroups.78-nup.15.png){ width=30% }
![](tree_phenotypes_regression.out.dLFE.by.taxgroups.78-nup.16.png){ width=30% }
![](tree_phenotypes_regression.out.dLFE.by.taxgroups.78-nup.17.png){ width=30% }
![](tree_phenotypes_regression.out.dLFE.by.taxgroups.78-nup.18.png){ width=30% }
![](tree_phenotypes_regression.out.dLFE.by.taxgroups.78-nup.19.png){ width=30% }
![](tree_phenotypes_regression.out.dLFE.by.taxgroups.78-nup.20.png){ width=30% }


Coefficient of determination ($R^2$) and regression direction (red - positive slope, blue, negative slope) for GLS regression between different traits and mean ΔLFE in regions relative to CDS start and end, for different taxonomic subgroups. Significant values (*p*-value < 0.01) are marked with white dots.

**A** Genomic-GC%


## Supplementary Figure S12. Dependence of ΔLFE profiles on temperature

![](correlation_between_series.temp.combined.png){ width=60% }  
ΔLFE
![](heatmap_profile_legend.png){ width=20% }

**A** Correlation between ΔLFE calculated using standard temperature ($37^{\circ}\text{C}$) and native temperature (see methods), at each position relative to CDS start, for species grouped by native temperature range. Correlations were calculated for a random sample (*N*=71) of species (bacteria and archaea) for which native temperature data is available.
**B** Comparison of individual mean ΔLFE profiles using calculated using standard temperature ($37^{\circ}\text{C}$) and native temperature.

