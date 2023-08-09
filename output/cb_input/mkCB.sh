#!/usr/bin/env bash

###CMD to run:
# bash mkCB.sh meta.csv


### MODIFY ###

#set to parent dir
export CBDATAROOT='/pl/active/dow_lab/dylan/eq_synovial_scRNA/analysis/output/cb_input'

outDir="../cb_output"
cmd1="cbBuild -r -i ./cellbrowser.conf -o ${outDir}"

#generate configs:

### setting:
name="equine_sfn3n3"
tags="10x"
organism="Equus caballus"
project="Equine synovial fluid atlas"
dotSize=2

title="A single-cell RNA sequencing atlas of equine synovial fluid in health and OA"
submitter="Dylan Ammons"
version=1
submission_date="2023-08-09"

### END MODIFY ###

while read line
do

    subName=( $(echo $line | cut -f1 -d',' --output-delimiter=' ') )
    defaultClus=( $(echo $line | cut -f2 -d',' --output-delimiter=' ') )
    pertyName=( $(echo $line | cut -f3 -d',' --output-delimiter=' ') )    

    echo 'shortLabel= '${pertyName} > $subName/cellbrowser.conf

    echo 'enumFields = ["'${defaultClus}'"]' >> $subName/cellbrowser.conf ### modify this to make universal

    echo 'coords=[{"file":"umap.coords.tsv", "shortLabel":"UMAP"}, ]' >> $subName/cellbrowser.conf


    echo 'clusterField = "'${defaultClus}'"' >> $subName/cellbrowser.conf
    echo 'labelField = "'${defaultClus}'"' >> $subName/cellbrowser.conf


    echo 'organisms = ["'${organism}'"]' >> $subName/cellbrowser.conf
    echo 'projects = ["'${project}'"]' >> $subName/cellbrowser.conf
    echo "radius = ${dotSize}" >> $subName/cellbrowser.conf

    echo 'quickGenesFile = "quickGenes.csv"' >> $subName/cellbrowser.conf

    #immutables:
    echo 'exprMatrix="exprMatrix.tsv.gz"' >> $subName/cellbrowser.conf
    echo 'geneIdType="raw"' >> $subName/cellbrowser.conf
    echo 'meta="meta.tsv"' >> $subName/cellbrowser.conf
    echo 'markers=[{"file":"markers.tsv", "shortLabel":"Cluster-specific markers"}]' >> $subName/cellbrowser.conf
    echo 'unit = "log of read count/UMI"' >> $subName/cellbrowser.conf
    echo 'matrixType="auto"' >> $subName/cellbrowser.conf
    
done < $1

#make main conf file
echo 'name = "'${name}'"' > cellbrowser.conf
echo 'shortLabel = "'${project}'"' >> cellbrowser.conf
echo 'tags = ["'${tag}'"]' >> cellbrowser.conf

#make main desc.conf
echo 'title = "'${title}'"' > desc.conf
echo 'submitter = "'${submitter}'"' >> desc.conf
echo 'version = "'${version}'"' >> desc.conf
echo 'submission_date = "'${submission_date}'"' >> desc.conf

#run cb
echo -e "\t$ ${cmd1}"
time eval $cmd1

cp $outDir/index.html $outDir/index.aspx