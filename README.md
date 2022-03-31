# Inference and Embedding of Multi-modal Networks

## Study information
Here, we present the code for our two-step framework for the analysis of multi-modal data, capturing the interplay between the different information layers. This framework consisted of inferring a multi-modal network and embedding the nodes into a low dimensional space for the effective exploration of similarities between nodes and data modalities.

<img src="https://user-images.githubusercontent.com/50441257/160692229-d5c2366e-880d-42a0-a910-30effc067aef.png" width="700" >


Code for the manuscript **Network Embedding across Multiple Tissues Elucidates Multi-modal Context of Host Factors Important for COVID-19 Infection' Yue Hu,  Ghalia Rehawi, Lambert Moyon, Christoph Ogris, Janine Knauer-Arloth, Florian Bittner, Annalisa Marsico,  Nikola S. Mueller (2022)**


## Data
Openly available public data Genotype-Tissue Expression (GTEx) from https://gtexportal.org/home/
complemented by confidential data on genotypes and phenotypes which cannot be disclosed here. Special access can be granted by application to NCBI dbGAP Portal (https://dbgap.ncbi.nlm.nih.gov/).

## Conda environments
1. Install anaconda and snakemake. (We used conda 4.11.0; snakemake-minimal=5.32.1=py_0)
2. Conda environment with all packages are to be created automatically by snakemake at each step 

## Code
Snakemake workflow manager was used for the different steps of the analysis:
1. Preprocessing
2. Polygenic Risk Score (PRS) calculation
3. KiMONo
4. Network embedding
5. Analysis & Figures

Same directory structure can be found for Snakemake files and code scripts.
Instructions on execution of snakemake workflows can be found in each directory in form of README.txt.
