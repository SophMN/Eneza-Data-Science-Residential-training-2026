#!/bin/bash

#for file in GSM*; do
#sample=${file%%_*}
#mkdir -p "$sample"
#done

#for file in GSM*.gz; do
#sample=${file%%_*}
#echo mv "$file" "$sample/"
#done

for dir in GSM*/; do
cd "$dir" || exit
mv *_barcodes.tsv.gz barcodes.tsv.gz
mv *_features.tsv.gz features.tsv.gz
mv *_matrix.mtx.gz matrix.mtx.gz
cd ..
done
