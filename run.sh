#!/usr/bin/python
inputDir=dev
outputDir=devout
mkdir -p $inputDir $outputDir
python fragment_forest.py --data_dir $inputDir --sub_id 0 --nodes node82 --save_dir $outputDir --stop stop_words.txt --lemma lemma.txt
