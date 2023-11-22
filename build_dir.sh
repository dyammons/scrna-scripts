#!/bin/bash


mkdir -p ../input/ ../output/ ../output/s1/ ../output/s2/  ../output/s3/  ../output/clustree/  ../output/singleR/  ../output/cb_input/  ../output/cb_output/  ../output/viln


for dir in $@
do
  baseDir=../output/$dir
  mkdir -p $baseDir $baseDir/pseudoBulk/ $baseDir/linDEG ../output/viln/$dir/
done
