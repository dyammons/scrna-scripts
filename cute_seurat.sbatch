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
