#!/usr/bin/env bash

###CMD to run:
# bash mkCB.sh ./metaData/meta.csv

#set to parent dir
export CBDATAROOT='../output/cb_input/'

outDir="../output/cb_output/"
cmd1="cbBuild -i ./cellbrowser.conf -o ${outDir}" #use this if you have more than 1 subset
# cmd1="cbBuild -r -i ./cellbrowser.conf -o ${outDir}" #use this if you have more than 1 subset

#generate configs:

### setting:
name="expt_1"
tags="10x"
organism="Canis lupus familiaris"
project="Single cell atlas"
dotSize=3

title="A single-cell RNA sequencing atlas of ...."
submitter="Dylan Ammons"
version=1
submission_date="2023-07-25"


while read line
do

    subName=( $(echo $line | cut -f1 -d',' --output-delimiter=' ') )
    defaultClus=( $(echo $line | cut -f2 -d',' --output-delimiter=' ') )
    pertyName=( $(echo $line | cut -f3 -d',' --output-delimiter=' ') )    

    echo 'shortLabel= '${pertyName} > ../output/cb_input/$subName/cellbrowser.conf

    echo 'enumFields = ["'${defaultClus}'"]' >> ../output/cb_input/$subName/cellbrowser.conf ### modify this to make universal

    echo 'coords=[{"file":"umap.coords.tsv", "shortLabel":"UMAP"}, ]' >> ../output/cb_input/$subName/cellbrowser.conf


    echo 'clusterField = "'${defaultClus}'"' >> ../output/cb_input/$subName/cellbrowser.conf
    echo 'labelField = "'${defaultClus}'"' >> ../output/cb_input/$subName/cellbrowser.conf


    echo 'organisms = ["'${organism}'"]' >> ../output/cb_input/$subName/cellbrowser.conf
    echo 'projects = ["'${project}'"]' >> ../output/cb_input/$subName/cellbrowser.conf
    echo "radius = ${dotSize}" >> ../output/cb_input/$subName/cellbrowser.conf

    echo 'quickGenesFile = "quickGenes.csv"' >> ../output/cb_input/$subName/cellbrowser.conf

    #immutables:
    echo 'exprMatrix="exprMatrix.tsv.gz"' >> ../output/cb_input/$subName/cellbrowser.conf
    echo 'geneIdType="raw"' >> ../output/cb_input/$subName/cellbrowser.conf
    echo 'meta="meta.tsv"' >> ../output/cb_input/$subName/cellbrowser.conf
    echo 'markers=[{"file":"markers.tsv", "shortLabel":"Cluster-specific markers"}]' >> ../output/cb_input/$subName/cellbrowser.conf
    echo 'unit = "log of read count/UMI"' >> ../output/cb_input/$subName/cellbrowser.conf
    echo 'matrixType="auto"' >> ../output/cb_input/$subName/cellbrowser.conf
    
done < $1

#make main conf file
echo 'name = "'${name}'"' > ../output/cb_input/cellbrowser.conf
echo 'shortLabel = "'${project}'"' >> ../output/cb_input/cellbrowser.conf
echo 'tags = ["'${tag}'"]' >> ../output/cb_input/cellbrowser.conf

#make main desc.conf
echo 'title = "'${title}'"' > ../output/cb_input/desc.conf
echo 'submitter = "'${submitter}'"' >> ../output/cb_input/desc.conf
echo 'version = "'${version}'"' >> ../output/cb_input/desc.conf
echo 'submission_date = "'${submission_date}'"' >> ../output/cb_input/desc.conf

#run cb
echo -e "\t$ ${cmd1}"
time eval $cmd1

cp $outDir/index.html $outDir/index.aspx