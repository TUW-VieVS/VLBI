#!/bin/bash

declare -a folders=("DATA/NGS" "DATA/vgosDB" "DATA/MASTER" "TRP/VMF1" "TRP/VMF3" "TRP/GRAD" "TRP/LHG" "TRP/RAYTRACING_DATA" "ATM/APL_VIENNA" "ATM/GSFC")

for i in "${folders[@]}" 
do
	echo "link $i"
	rm -r ../../../VLBI/$i/*  > /dev/null 2>&1
	ln -sfn /home/members/vievs/VLBI/$i/* ../../../VLBI/$i
done


declare -a files=("EOP")
for i in "${files[@]}" 
do
	echo "link $i"
	rm -r ../../../VLBI/$i/*.txt  > /dev/null 2>&1
	ln -sfn /home/members/vievs/VLBI/$i/*.txt ../../../VLBI/$i
done






