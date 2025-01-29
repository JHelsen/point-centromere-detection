# PCAn v0.1

Point Centromere Annoation (PCAn) v0.1 is a tool used to identify and annotate point centromeres from the Saccharomycetaceae clade. It is designed in Python to work on genome assemblies at any level of completeness (contig, scaffold, chromosome, or complete).

## Installation
This tool works in the following systems

<img src="https://github.com/primefaces/primeicons/blob/master/raw-svg/apple.svg" width="50" height="50"> <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/35/Tux.svg/1280px-Tux.svg.png" width="40" height="40">

Miniconda was used as the environment and package manager for this protject. To set up Miniconda, please refer to the [documentation](https://docs.anaconda.com/miniconda/index.html)
  
  1. Download the PCAn.zip file using this [link](https://github.com/JHelsen/point-centromere-detection/blob/main/PCAn/PCAn.zip).
  2. Navigate to the downloaded .zip file and extract all files into a separate folder
  3. Navigate to the extracted folder on your terminal and use the following command to create the environment
     ```bash
      conda env create -n pcan --file pcan_specs.yml
     ```
  4. Activate the environment using the following command prior to using the tool
     ```bash
     conda activate pcan
     ``` 

## Using PCAn
The PCAn zip archive comes with all the files required to run PCAn. PCAn uses [MEME](https://meme-suite.org/meme/) motifs that have been calculated for each clade within the Saccharomycetaceae to identify and extract centromeres from a given genome assembly. If the user allows it, PCAn can also use a combination of BLAST and a [Python version of NCBI ORFFinder](https://github.com/Chokyotager/ORFFinder) to identify genes flanking the identified centromere, and perform a synteny check against the ancestral gene order (calculated and reported in earlier publications [1](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000485)) to further confirm the veracity of the extracted centromeres. 

PCAn requires the following user inputs 
  1. The path to the genome assembly
  2. Selection of the motif that is phylogenetically closest to the input genome assembly
  3. Permission to perform the synteny check using BLAST ("Yes" or "No") - this adds to the computing time

## PCAn outputs
You can find sample predictions for _Saccharomyces cerevisiae_ in the "sample_outputs" folder

### CENsequences.txt
This text file contains the table of centromere sequences identified and extracted by PCAn. The table also contains FIMO scores and overall scores. The details on how the overall scores are calculated can be found in our preprint.
It also comes with other information below the table after some *s - number of putative centromeres, path to the genome used, motifs and thresholds used  

### SyntenyPlot.png
The first number corresponds to the ancestral chromosome, L or R refer to the ancestral chromosome arm, and the last number indicates how far the protein lies from the ancestral centromere (with 1 right next to the ancestral centromere). Proteins on left arms are colored in blue, and proteins on right arm are colored in pink. Contigs are indicated on the right.

### Please note
Motifs and thresholds are optimized to ensure low false positive and negative rates, and to optimize computing time. Most of the false negatives can be found by lowering the CDEIII detection threshold. However, this will increase computing time and the number of false positive hits.



## Please cite

>J. Helsen, K. Ramachandran, G. Sherlock, G. Dey, Centromeres evolve progressively through selection at the kinetochore interface. bioRxiv [Preprint] (2025). https://doi.org/10.1101/2025.01.16.633479.



## Contact
If you run into any issues while running PCAn, please either raise an issue here or contact Jana Helsen
