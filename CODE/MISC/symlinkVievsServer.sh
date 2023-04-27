#!/bin/bash

declare -a folders=("DATA/NGS" "DATA/vgosDB" "DATA/MASTER" "TRP/VMF1" "TRP/VMF3" "TRP/GRAD" "TRP/RAYTRACING_DATA" "ATM/APL_VIENNA" "NTSL/NTAL/VIE")

for i in "${folders[@]}" 
do
	echo "link $i"
	rm -r ../../$i/*  > /dev/null 2>&1
	ln -sfn /home/members/VieVS/VLBI/$i/* ../../$i
done


declare -a files=("EOP")
for i in "${files[@]}" 
do
	echo "link $i"
	rm -r ../../$i/*.txt  > /dev/null 2>&1
	ln -sfn /home/members/VieVS/VLBI/$i/*.txt ../../$i
done






