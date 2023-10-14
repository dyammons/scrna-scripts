# Single-cell RNA sequencing skeleton for mulit-sample, multi-condition R analysis
This is a directory skeleton with code templates for scrna-seq analysis designed to process data for mulit-sample, multi-condition experiments. Useage is intended to facilitate the generation of uniform preliminary analysis across projects and promote an understanding of the computational steps involved in scRNA-seq analysis.  

To use this repository, you will need to clone the repo to your working space, create a computing environement with required software, and provide count matricies as an input.

There are several ways to clone the repository (including through GitHub Desktop), but here is the command to do it through a linux terminal. If working on Alpine, it would be ideal to clone this to a `Peta Library allocation`, but if that is not avalible `Scratch` or `Projects` directories will also work. If using `Scratch` please note that data is deleted every 90 days, so you will need to complete regular backups and if using `Projects` there is limited storage, so there may not be sufficent room to complete the analysis. 

```sh
git clone https://github.com/dyammons/scrna-seq-skeleton-1.git
```

For standardization of analysis it is reccomended that all analysis be complete on the Alpine supercomputer in a custom Conda environment; however, all required software and code can be installed/run locally in R Studio if prefered.
To set up the computing enviornment a `.yml` file containing all the information required to build a functional Conda env is provided in the `scRNA_r_env.yml` file.  
In the cloned repository, run the follinw code to generate the env:

```sh
conda env create -n scRNA_r_env -f scRNA_r_env.yml ### need to add .yml and test functionality
```
Note: Installation of software can be completed manually using Conda and steps to create a suitable environment are avalible [here](). (TO DO: add actual link....)

For ease of creating required output directories the `build_dir.sh` file is provided. This short script will generate the expected output directory with subfolders for each major cell type subset.  
Run with:
```sh
bash build_dir.sh allCells tcell bcell #where each argument is the name of cell subtype - change the example as needed
```

The final step to get started is to put your single-cell count matrices in the `input` directory. File structure within `input` should be that each sample has its own directory with the corresponding `features.tsv.gz`, `matrix.mtx.gz`, and `barcodes.tsv.gz`.

From there you should be able to load the data in with `script1.R` and process some data!

Questions? Submit an issue or reach out to Dylan Ammons directly.
