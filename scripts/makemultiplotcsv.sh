#!/bin/bash
mkdir recoverymultiplot2
touch recoverymultiplot2/filenames
mkdir recoverymultiplot2/rfamonalien
mkdir recoverymultiplot2/alienonrfam
echo "" > recoverymultiplot2/filenames
counter=1
for i in structuredalienseedoutput4-*.csv; do 
	echo "recoverymultiplot2/rfamonalien/$counter.sorted.csv " >> recoverymultiplot2/filenames
	echo "$i"
	cut -d $'\t' -f 13 $i > recoverymultiplot2/rfamonalien/$counter.csv; 
	cut -d $'\t' -f 13 $i > recoverymultiplot2/alienonrfam/$counter.csv;
	sort -k 1 -n recoverymultiplot2/rfamonalien/$counter.csv > recoverymultiplot2/rfamonalien/$counter.sorted.csv;
	sort -k 1 -n recoverymultiplot2/alienonrfam/$counter.csv > recoverymultiplot2/alienonrfam/$counter.sorted.csv;
	counter=$[$counter +1]
done
usedfilenames=$(<recoverymultiplot2/filenames)
echo "usedfilenames"
pr -mts $usedfilenames > recoverymultiplot2/rfamonalien/allsorted
pr -mts $usedfilenames > recoverymultiplot2/alienonrfam/allsorted

