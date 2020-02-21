Perform analysis of mRNA folding bias (Local Folding Energy) based on randomization experiments.
See publication (Peeri and Tuller 2020) for description of methods.

Important source files:

* taxonomy.py -- Create multiple-species plots, based on taxonomic tree or PCA

* profile_plots.py -- Create plot for individual species

* tree_traits_effects_analysis_with_taxgroups.r -- Perform GLS regression and other analyses for multiple species based on phylogenetic tree

* config.py -- Configuration parameters (including database connections, etc.)

* r-rnafold -- Singularity recipe for initializing container with R environment

