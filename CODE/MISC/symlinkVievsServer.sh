#!/bin/bash

declare -a folders=("DATA/NGS" "DATA/vgosDB" "DATA/MASTER" "TRP/VMF1" "TRP/VMF3" "TRP/GRAD" "TRP/RAYTRACING_DATA" "NTSL/NTAL/VIE")

for i in "${folders[@]}" 
do
	echo "link $i"
	rm -r ../../$i/*  > /dev/null 2>&1
	ln -sfn /data/DATA4VLBI/$i/* ../../$i
done


declare -a files=("EOP")
for i in "${files[@]}" 
do
	echo "link $i"
	rm -r ../../$i/*.txt  > /dev/null 2>&1
	ln -sfn /data/DATA4VLBI/$i/*.txt ../../$i
done






