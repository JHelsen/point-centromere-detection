The scripts in this directory describe how to setup the conda environment and use the Python scripts in the "scripts" directory to reconstruct the species tree from our study.

## Setting up the Miniconda environment
Miniconda was used as the environment and package manager for this protject. To set up Miniconda, please refer to the [documentation](https://docs.anaconda.com/miniconda/index.html)

You can setup the conda environment using the following command. You would need to download the files in this directory prior to this step.
```bash
conda env create -f /path/to/homology_searches_specs.yml
```

Before running the scripts, activate the environment using
```bash
conda activate homology_searches
```

## Selecting marker sequences 
We selected 1403 markers from [Opulente et al. 2024](https://doi.org/10.1126/science.adj4503). The relevant files from this study can be found [here](https://plus.figshare.com/articles/dataset/Dataset_supporting_Figure_2_Phylogenomics/22806416?backTo=/collections/Genomic_and_ecological_factors_shaping_specialism_and_generalism_across_an_entire_subphylum/6714042). 

Prior to using the marker sequences, we used the scripts in "scripts/marker-selection" to identify and exclude any marker sequence that
1. did not have a homolog in _S. cerevisiae_ - 113 markers
2. redundantly mapped with another marker to the same _S. cerevisiae_ protein - 13 markers

The scripts were used in the following order. The relevant documentation for each script can be found in the comments within the script itself:
1. phmmer_run.py - to use phmmer to map each marker sequence back to a protein in _S. cerevisiae_
2. phmmer_extract_top_hit.py - to extract the top hit from each phmmer search, and identify markers which mapped back to the same protein

## Retrieving homologs for the selected markers
Prior to this step, all genomes were organized into a single folder. BLAST databases for these genomes were built using 
```bash
makeblastdb -in <genome_name>.fasta -out <genome_name> -dbtype nucl
```
These BLAST databases were then moved to a separate directory.

The relevant documentation for each script can be found in the comments within the script itself.
1. [1_parallelized_tblastn_searches.py](https://github.com/JHelsen/point-centromere-detection/blob/main/species-tree-reconstruction/scripts/1_parallelized_tblastn_searches.py) - to perform tBLASTn searches for each selected marker against the constructed BLAST databases
2. [2_tblastn_results_analysis_orf_extraction.py](https://github.com/JHelsen/point-centromere-detection/blob/main/species-tree-reconstruction/scripts/2_tblastn_results_analysis_orf_extraction.py) - to analyse the tBLASTn search results, reject all hits whose E-value is lower than the threshold of 1E-10, extract scaffolds, run ORFFinder, and write out the predicted ORFs
3. [3_reverse_blastp_against_cerevisiae.py](https://github.com/JHelsen/point-centromere-detection/blob/main/species-tree-reconstruction/scripts/3_reverse_blastp_against_cerevisiae.py) - to perform BLASTp of each selected ORF against the _S. cerevisiae_ proteome
4. [4_reverse_blast_analysis_for_species_tree.py](https://github.com/JHelsen/point-centromere-detection/blob/main/species-tree-reconstruction/scripts/4_reverse_blast_analysis_for_species_tree.py) - to analyse the BLASTp results, select those ORFs which correctly mapped back to the corresponding protein in the _S. cerevisiae_ proteome, and write those ORFs out separately

## Aligning the marker homolog datasets, building gene trees, and removing outliers
At this stage, you will have a folder created whose name ends with "_seqs_for_tree_building". This folder will contain subfolders, one for each input marker sequence. Each folder will contain a single FASTA file containing all the identified homologs for that marker.

5. [5_align_trim_tree_building.py](https://github.com/JHelsen/point-centromere-detection/blob/main/species-tree-reconstruction/scripts/5_align_trim_tree_building.py) - run this script inside each subfolder within the "_seqs_for_tree_building" folder to build alignments, trim them, and construct maximum-likelihood trees
6. [6_identifying_outliers_tree_markers.py](https://github.com/JHelsen/point-centromere-detection/blob/main/species-tree-reconstruction/scripts/6_identifying_outliers_tree_markers.py) - run this script within the "_seqs_for_tree_building" folder to systematically identify outliers - sequences whose branch lengths in the constructed tree is greater than 20 times the median branch length - within each subfolder. The code will also identify homolog sets which cannot automatically be handled for manual checking

Repeat 5. and 6. till no further outliers can be identified. Following this,
## Building species trees
7. [7_alignment_modification_for_species_tree_building.py](https://github.com/JHelsen/point-centromere-detection/blob/main/species-tree-reconstruction/scripts/7_alignment_modification_for_species_tree_building.py) - run this script within the "_seqs_for_tree_building" folder to import all trimmed alignments, rename each sequence within each alignment, and write the renamed alignments into a separate folder

Within this folder, use Goalign to concatenate the alignments
```bash
goalign concat *.phy > supermatrix_alignment.fasta
```

The concatenated supermatrix alignment is used to build the species tree using IQ-TREE
```bash
iqtree2 --seqtype AA -B 1000 -alrt 1000 --boot-trees --wbtl -m LG+G4 -mwopt --threads-max 24 -T AUTO -s <supermatrix_alignment.fasta>
```
## Dating the constructed species tree
8. [8_random_subsampling_of_supermatrix_alignment_sites.py](https://github.com/JHelsen/point-centromere-detection/blob/main/species-tree-reconstruction/scripts/8_random_subsampling_of_supermatrix_alignment_sites.py) - To make the tree dating more tractable, the supermatrix alignment was randomly subsampled and 3 sets of 10000 sites were extracted.

Each subsampled alignment was used to build another species tree using IQ-TREE. The topology was fixed using the species tree obtained earlier.
```bash
iqtree2 --seqtype AA -B 1000 --boot-trees --wbtl -m LG+G4 --score-diff ALL --threads-max 32 -T AUTO -mwopt -s <subsample_1_sequence.fasta> -g <fixed_topology.tree> -pre <fixed_topology.tree>
```
The 3 resulting trees were dated using MEGA as described in our study.
