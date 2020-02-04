#!/bin/bash
echo "user name:"
read user

echo "years: (e.g.: \"2012\" or \"*\" if you want all years)"
read years

echo "download NGS files? (y for yes)"
read ngs

echo "download vgosDB files? (y for yes)"
read vgosDB

echo "download RATRACING_DATA files? (y for yes)"
read radiate

echo "download EOPs"
scp $user@vievs:/home/members/VieVS/VLBI/EOP/*.txt ../../EOP/
echo "download APL"
scp $user@vievs:/home/members/VieVS/VLBI/ATM/APL_VIENNA/y$years.apl_r ../../ATM/APL_VIENNA
if [ "$ngs" = "y" ]; then
	echo "download NGS"
	rsync -avzh $user@vievs:/home/members/VieVS/VLBI/DATA/NGS/$years* ../../DATA/NGS
else
	echo "skipping NGS"
fi
if [ "$vgosDB" = "y" ]; then
	echo "download vgosDB"
	rsync -avzh $user@vievs:/home/members/VieVS/VLBI/DATA/vgosDB/$years* ../../DATA/vgosDB
else
	echo "skipping vgosDB"
fi
if [ "$radiate" = "y" ]; then
        echo "download RAYTRACING_DATA"
        rsync -avzh $user@vievs:/home/members/VieVS/VLBI/TRP/RAYTRACING_DATA/$years* ../../TRP/RAYTRACING_DATA
else
        echo "skipping RAYTRACING_DATA"
fi
echo "download VMF1"
scp  $user@vievs:/home/members/VieVS/VLBI/TRP/VMF1/y$years.vmf1_r ../../TRP/VMF1
echo "download VMF3"
scp  $user@vievs:/home/members/VieVS/VLBI/TRP/VMF3/y$years.vmf3_r ../../TRP/VMF3
echo "download GRAD"
scp  $user@vievs:/home/members/VieVS/VLBI/TRP/GRAD/y$years.grad_r ../../TRP/GRAD
