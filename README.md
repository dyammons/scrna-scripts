# Single-cell RNA sequencing repository for preliminary R analysis
This repository contains the instructions to create a directory structure for analysis and pull down Singularity container to run a Seurat pipeline to integrate multiple samples into one dataset.
Useage is intended to facilitate the generation of uniform preliminary analysis across projects and promote an understanding of the computational steps involved in scRNA-seq analysis.  

To use this repository, you will need to clone the repo to your working space, pull down a Singularity container a computing environement with required software, and provide count matricies as an input.

<br>

## Set up the directory structure and get scripts
There are several ways to clone the repository (including through GitHub Desktop), but here is the command to do it through a linux terminal. If working on Alpine, it would be ideal to clone this to a `Peta Library allocation`, but if that is not avalible `Scratch` or `Projects` directories will also work. If using `Scratch` please note that data is deleted every 90 days, so you will need to complete regular backups and if using `Projects` there is limited storage, so there may not be sufficent room to complete the analysis. 

```sh
mkdir -p /scratch/alpine/$USER/scrna-analysis/
cd /scratch/alpine/$USER/scrna-analysis/

git clone https://github.com/dyammons/scrna-scripts.git
```

Navigate into the repo:
```sh
cd scrna-scripts
```

For ease of creating the required output directories the `build_dir.sh` file is provided. This short script will generate the necessary output directories and subfolders for each major cell type.  
Run with:
```sh
bash build_dir.sh allCells tcell bcell #where each argument is the name of cell subtype - change the example as needed - you can always add more later by rerunning this script
```

Go up a level and you should now see `input` and `output` in addition to the original `scrna-scripts` directory.
```sh
cd ..
ls

```

The directory stucture in `output` will look something like this:
```sh
output/
├── allCells
│   ├── linDEG
│   └── pseudoBulk
├── cb_input
├── cb_output
├── clustree
├── linDEG
├── s1
├── s2
├── s3
├── singleR
└── viln
    └── allCells
```

## Bring count matricies into the input directory

You will now need to copy your single-cell count matrices in the `input` directory. File structure within `input` should be such that each sample has its own directory with the corresponding `features.tsv.gz`, `matrix.mtx.gz`, and `barcodes.tsv.gz` (dir tree below).
```sh
input/
├── sample1
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── sample2
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── sample3
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── ...
```

It is best to avoid moving files manually, so here is the approach I commonly use.
```sh
#set string array with names of dirs you want to get data from
dirs=$( ls -l /scratch/alpine/$USER/project_01/02_scripts | grep "^d" | awk '{print $9}' )
declare -a StringArray=($dirs)
```

```sh
cd input

#loop through the array create sample sub-directories thene copy the filtered_feature_bc_matrix
for val in "${StringArray[@]}"; do
  folder="./input/$val/"
  mkdir $folder
  
  filez="/scratch/alpine/dyammons@colostate.edu/proj03_k9_duod/02_scripts/$val/outs/filtered_feature_bc_matrix/*"
  cp $filez $folder
done
```

<details><summary>Code as a script</summary>
<p>

```sh
#!/usr/bin/env bash

###MODIFY as needed!
###Useage: bash getData.sh
###Run this in the main analysis directory (or change the paths in the code as needed.

#set string array with names of dirs you want to get data from
dirs=$( ls -l /scratch/alpine/$USER/project_01/02_scripts | grep "^d" | awk '{print $9}' )
declare -a StringArray=($dirs)
#mkdir input

#loop through the array create sample sub-directories thene copy the filtered_feature_bc_matrix
for val in "${StringArray[@]}"; do
  folder="./input/$val/"
  mkdir $folder
  
  filez="/scratch/alpine/dyammons@colostate.edu/proj03_k9_duod/02_scripts/$val/outs/filtered_feature_bc_matrix/*"
  cp $filez $folder
done
```

</p>
</details>

<br> 

## Collect the software container

With the input datasets in place we are nearly ready to get the code running! The next step is to get the `Singularity` container we will be using to run the script.  
So, let's pull it down from syslabs:
```sh
#establish the connection to syslabs
apptainer remote add --no-login SylabsCloud cloud.sycloud.io
apptainer remote use SylabsCloud
```
```sh
#move into the scripts dir and pull down the sif
cd scrna-scripts/
singularity pull --arch amd64 library://dyammons/r-env/r4.3.1-seurat:v1
```

<details><summary>If the above fails you can cp a copy from my scratch space.</summary>
<p>
  
```sh
cd scrna-scripts/
cp /scratch/alpine/dyammons@colostate.edu/dump/scrna-scripts/r4.3.1-seurat_v1.sif .
```

</p>
</details>

Let's make sure we can enter the container and the software is accessible for our use.  
To do this we will launch a shell. This is very similar to what `conda activate env` would do if you are familiar with `conda`.
```sh
#it is important to bind (-B) a directory at least 1 level up from the scripts folder
singularity shell -B $PWD/../ r4.3.1-seurat_v1.sif
```

While in the container we have access to all the software. So, let's launch an R session to ensure we can `source` the `custonFunctions.R` file that will be key to running the code.
```sh
R
source("./customFunctions.R")
```

If all the packages load in no problem, then we are good to move forward!

Since we are already in the container, let's run the code to generate the QC paramters so we can set threshold for the pipeline.
```r
load10x(din = "../input/", dout = "../output/s1/", outName = "qc_test", testQC = T)
#Saving 7 x 7 in image
```

We can now use our file navigator panel to inspect the QC plots (`../output/s1`).

View the files and decide on thresholds. Err on the side of caution and set them permissively as we can always go back and increase the stringency later on.

We will code in the thresholds by opening the `script1.R` file and customizing the `MODIFY` section of the script.
```r
######### MODIFY #########

#set output name -- recommend including data and sample size
experiment <- "pbmc_analysis_110723"
outName <- "allCells"
testQC <- F
nFeature_RNA_high <- 5500
nFeature_RNA_low <- 100
percent.mt_high <- 12.5
nCount_RNA_high <- 60000
nCount_RNA_low <- 200

########## END MODIFY #########
```

With the sif downloaded we are ready to run the code. Before submitting a job,

Lastly, we will enter in some metadata that will be used to colorize the samples and get short sample names loaded in.

To do this we will open the `./metaData/refColz.csv` in a text editor and modify it as desired.
The columns `orig.ident` should exactly match the samples names as defined in the `input` sub-directories. The `name` column can be anything you want, typically a short hand for the sample name.

Once the values are entered in the Rscript and the metadata is entered we are ready to run the preliminary script.

```sh
#get out of the container
q()
n
exit
```

Run a job with this script. Run from inside the scripts directory.
```sh
#!/usr/bin/env bash

#SBATCH --job-name=seu_prelim
#SBATCH --ntasks=10       # 10 worked well for  6 samples with ~5k cells each, scale up if more samples
#SBATCH --nodes=1         # this script is designed to run on one node
#SBATCH --time=02:00:00   # set time; default = 4 hours

#SBATCH --partition=amilan  # modify this to reflect which queue you want to use. Either 'shas' or 'shas-testing'
#SBATCH --qos=normal      # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'

#SBATCH --mail-type=END   # Keep these two lines of code if you want an e-mail sent to you when it is complete.
#SBATCH --mail-user=dyammons@colostate.edu ### change to your email ###

#SBATCH --output=seu_prelim_%j.log  #modify as desired - will output a log file where the "%j" inserts the job ID number

######### Instructions ###########

#run R script
singularity exec -B  r4.3.1-seurat.sif Rscript script1.R

```
Questions? Submit an issue or reach out to Dylan Ammons directly.
