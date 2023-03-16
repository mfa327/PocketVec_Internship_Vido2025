# PocketVec

![alt text](https://gitlabsbnb.irbbarcelona.org/acomajuncosa/pocketvec/blob/master/data/png_images/PocketVec.png?raw=true)

Pocket descriptors characterize protein binding sites in the shape of numerical vectors, which enable the high-throughput exploration of the pocket space.

**PocketVec** is a novel strategy to generate protein binding site descriptors based on the assumption that similar pockets bind similar ligands. The approach is built upon inverse virtual screening: the prioritization of a predefined set of small molecules is expected to be more correlated between a pair of similar pockets than between a pair of dissimilar ones. In short, each pocket is represented by the set of numerical ranks provided by small molecules in the docking calculations. Further conceptual and methodological details are best described in the original PocketVec publication: REF. 


The **PocketVec Repository** holds the code needed to create a PocketVec descriptor for any protein binding site of interest together with all results presented along the manuscript. The directory structure is specified as follows:

* `code`: funcions, scripts and notebooks to generate PocketVec descriptors. 
* `data`: sets of small molecules to dock. (?)
* `examples`: example exercises. Check section `Basic usage` for further details. 
* `results`: PocketVec paper results.



# Installation

1. Clone this repository to your local PocketVec folder:
        
        cd ~ && mkdir -p pocketvec && cd pocketvec
        git clone https://gitlabsbnb.irbbarcelona.org/acomajuncosa/pocketvec.git

2. Create a conda environment with all the requirements:

        conda env create --name pocketvec_env --file=environment.yml
        conda activate pocketvec_env


# Running PocketVec



# Speed



# Basic usage



# External dependencies


# Citation
