# Repository containing code and instructions to complete preliminary analysis of single-cell RNA sequencing data in R 
This repository contains the instructions to create the expected directory structure, pull down a Singularity container, and to run a Seurat pipeline that will integrate multiple samples into one dataset.
Purpose: facilitate the generation of uniform preliminary analysis across projects and promote an understanding of the computational steps involved in scRNA-seq analysis.  


## Order of operations:
To use this repository, you will need to clone the repo to your working space, pull down a Singularity container with required software pre-installed, and provide count matrices as an input.

1. [Get the scripts](#set-up-the-directory-structure-and-get-scripts)
2. [Collect count matrices](#bring-count-matrices-into-the-input-directory)
3. [Collect the container](#collect-the-software-container)
4. [Test the software container](#test-the-software-container)
5. [Load in data and plot QC parameters](#load-in-data-and-plot-QC-parameters)
6. [Set thresholds and create a metadata file](#set-thresholds-and-create-a-metadata-file)
7. [Modify the provided `.sbatch` file and submit the job](#modify-the-provided-cute_seuratsbatch-file-and-submit-the-job)

<br>


## Set up the directory structure and get scripts
There are several ways to clone the repository we will be using today (including through GitHub Desktop), but here is the command to do it through a linux terminal.  

A general note: If working on Alpine, it would be ideal to clone this to a `Peta Library allocation`, but if that is not available `Scratch` or `Projects` directories will also work. If using `Scratch` please note that data is deleted every 90 days, so you will need to complete regular backups and if using `Projects` there is limited storage, so there may not be sufficient room to complete the analysis.  

For the purposes of today we will be working in `scratch` in a directory called `scrna-analysis`.

<br>

```sh
#make directory and navigate there
mkdir -p /scratch/alpine/$USER/scrna-analysis/
cd /scratch/alpine/$USER/scrna-analysis/
```

<br>

```sh
#clone the repo
git clone https://github.com/dyammons/scrna-scripts.git
```

<br>

```sh
#navigate into the repo
cd scrna-scripts
```

<br>

For ease of creating the required output directories the [`build_dir.sh`](./build_dir.sh) file is provided. This short script will generate the necessary output directories and subfolders for each major cell type.  

For this script you can enter multiple arguments where each argument is the name of cell subtype. Note: you can always add more later by rerunning this script

<br>

```sh
#run to create dir structure for "allCells"
bash build_dir.sh allCells 

#example with more cell types
#bash build_dir.sh allCells tcells bcells
```

<br>

Go up a level and you should now see `input` and `output` in addition to the original `scrna-scripts` directory.

<br>

```sh
cd ..

ls
#input  output  scrna-scripts
```

The directory structure in `output` will look something like this:
<details><summary>Show directory tree</summary>

<br>

```sh
output/
├── allCells
│   ├── linDEG
│   └── pseudoBulk
├── cb_input
├── cb_output
├── clustree
├── s1
├── s2
├── s3
├── singleR
└── viln
    └── allCells
```

</details>


## Bring count matrices into the input directory

You will now need to copy your single-cell count matrices into the `input` directory. File structure within `input` should be such that each sample has its own directory with the corresponding `features.tsv.gz`, `matrix.mtx.gz`, and `barcodes.tsv.gz` (dir tree below).

<details><summary>Show expected directory structure</summary>

<br>
    
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

</details>

<br> 

#### It is best to avoid moving files manually, so here is the approach I commonly use.

<br>

```sh
#navigate to input
cd input
```
<br>

Create a string array that contains the sample names.

<br>

```sh
#indicate path to directory containing the output files
path=/scratch/alpine/$USER/project_01/02_scripts

#set string array with names of dirs you want to get data from
dirs=$( ls -l $path | grep "^d" | awk '{print $9}' )
declare -a StringArray=($dirs)
```
<br>

Copy the data over.

<br>

```sh
#loop through the array to create sample sub-directories then copy the contents of filtered_feature_bc_matrix
for val in "${StringArray[@]}"; do
  folder="./$val/"
  mkdir $folder
  
  filez="$path/$val/outs/filtered_feature_bc_matrix/*"
  cp $filez $folder
done
```

<br>

<details><summary>Code as a script</summary>

<br>

Create a script file.
```sh
nano getData.sh
```

Copy the contents below then MODIFY paths as needed for your directory structure.
```sh

#!/usr/bin/env bash

###MODIFY as needed!
###Usage: bash getData.sh
###Run this in the input directory (or change the paths in the code as needed).


### User input ###

#indicate path to directory containing the output files
path=/scratch/alpine/$USER/project_01/02_scripts

### END User input ###



### CODE ###

#set string array with names of dirs you want to get data from
dirs=$( ls -l $path | grep "^d" | awk '{print $9}' )
declare -a StringArray=($dirs)

#loop through the array to create sample sub-directories then copy the filtered_feature_bc_matrix
for val in "${StringArray[@]}"; do
  folder="./$val/"
  mkdir $folder
  
  filez="$path/$val/outs/filtered_feature_bc_matrix/*"
  cp $filez $folder
done

### END CODE ###

```

Run the script in `input` to move the files to the required location.  

```sh
bash getData.sh
```

</details>

<br> 


## Collect the software container

With the input in place we are nearly ready to get the code running! The next step is to get the `Singularity` container we will be using to run the script.  
So, let's pull it down from Syslabs.

<br>

```sh
#move into the scripts dir and pull down the sif
cd ../scrna-scripts/
singularity pull --arch amd64 library://dyammons/r-env/r4.3.1-seurat:v1
```

<details><summary>If pull fails, try running this then the above code again.</summary>

<br>

```sh
#establish the connection to syslabs
apptainer remote add --no-login SylabsCloud cloud.sycloud.io
apptainer remote use SylabsCloud
```

</details>

<details><summary>If all of the above fails you can cp a copy from my scratch space.</summary>

<br>

```sh
cd scrna-scripts/
cp /scratch/alpine/dyammons@colostate.edu/dump/scrna-scripts/r4.3.1-seurat_v1.sif .
```

</details>


## Test the software container

Let's make sure we can enter the container and that the software is accessible for our use.  

To do this we will launch a `shell` to enter the container. This is very similar to what `conda activate env` would do if you are familiar with `conda`.

<br>

```sh
#it is important to bind (-B) a directory at least 1 level up from the scripts folder
singularity shell -B $PWD/../ r4.3.1-seurat_v1.sif
```

<br>

While in the container we have access to all the software. So, let's launch an `R` session to ensure we can `source` the `customFunctions.R` file that will be key to running the code.

<br>

```sh
R
source("./customFunctions.R")
```

<br>


## Load in data and plot QC parameters

If all the packages load in no problem, then we are good to move forward!

Since we are already in the container, let's run the code to generate the QC parameters so we can set thresholds for the pipeline.

<br>

```r
load10x(din = "../input/", dout = "../output/s1/", outName = "qc_test", testQC = T)
#Saving 7 x 7 in image
```

<br>

We can now use our file navigator panel to inspect the QC plots (`../output/s1`).

View the files and decide on thresholds.  

I recommend to err on the side of caution and set them permissively as we can always go back and increase the stringency later on.


## Set thresholds and create a metadata file

We will code in the thresholds by opening the `script1.R` file and customizing the `MODIFY` section of the script.

<br>

Excerpt provided here.
```r
######### MODIFY #########

#set output name
experiment <- "pbmc_analysis_20231129"
outName <- "allCells"

#set QC thresholds
nFeature_RNA_high <- 5500
nFeature_RNA_low <- 100
percent.mt_high <- 12.5
nCount_RNA_high <- 60000
nCount_RNA_low <- 200

########## END MODIFY #########
```

<br>

Lastly, we will enter in some metadata that will be used to colorize the samples and get short sample names loaded in.

To do this we will open the `./metaData/refColz.csv` in a text editor and modify it as desired.  

The columns `orig.ident` should exactly match the samples names as defined in the `input` sub-directories. The `name` column can be anything you want, typically a short hand for the sample name.

Once the values are entered in the Rscript and the metadata is entered we are ready to run the preliminary script. So, let's exit the container and prepare the `.sbatch` file.

<br>

```r
#quit the R session
q()
n
```

```sh
#leave the container
exit
```

<br>


## Modify the provided `cute_seurat.sbatch` file and submit the job

Open `cute_seurat.sbatch` in a text editor and modify it as desired.
Key parts to modify are:
 * `ntasks` the current default it set to 10. This worked will for 6 samples, may need to scale up running more samples.
 * `time` 2 hours should be good, but if running > 10 samples, may want to increase
 * `mail-user` change this one to your email so I don't get a notification that you ran a job (unless you want me to know)

<details><summary>Show script</summary>

<br>

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
#remove any loaded software
module purge

#run R script
singularity exec -B  $PWD/../ r4.3.1-seurat_v1.sif Rscript script1.R

```

</details>

<br>

```sh
#submit the job
sbatch cute_seurat.sbatch
```

<br>

The job should be competed in 1-3 hours depending on the number of samples you are integrating.

Questions? Submit an issue or reach out to Dylan Ammons directly.
