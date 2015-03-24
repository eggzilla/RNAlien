#!/bin/bash
mkdir recoverymultiplot2
mkdir recoverymultiplot2/rfamonalien
mkdir recoverymultiplot2/alienonrfam
counter=1
for i in structuredalienseedoutput4-*.csv; do 
	cut -d $'\t' -f 13 $i > recoverymultiplot2/rfamonalien/$counter.csv; 
	cut -d $'\t' -f 13 $i > recoverymultiplot2/alienonrfam/$counter.csv;
	sort -k 1 -n recoverymultiplot2/rfamonalien/$counter.csv > recoverymultiplot2/rfamonalien/$counter.sorted.csv;
	sort -k 1 -n recoverymultiplot2/alienonrfam/$counter.csv > recoverymultiplot2/alienonrfam/$counter.sorted.csv;
	counter=$[$counter +1]
done
pr -mts recoverymultiplot2/rfamonalien/*.sorted.csv > recoverymultiplot2/rfamonalien/allsorted
pr -mts recoverymultiplot2/alienonrfam/*.sorted.csv > recoverymultiplot2/alienonrfam/allsorted

