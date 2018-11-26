#!/bin/bash
echo "user name:"
read user

echo "years: (e.g.: \"2012\" or \"*\" if you want all years - you can use regular expressions)"
read years

echo "download NGS files? (y for yes)"
read ngs

echo "download vgosDB files? (y for yes)"
read vgosDB

echo "download EOPs"
scp $user@vievs:/home/members/vievs/VLBI/EOP/*.txt ../../../VLBI/EOP/
echo "download APL"
scp $user@vievs:/home/members/vievs/VLBI/ATM/APL_VIENNA/$years*.mat ../../../VLBI/ATM/APL_VIENNA
if [ "$ngs" = "y" ]; then
	echo "download NGS"
	rsync -avzh $user@vievs:/home/members/vievs/VLBI/DATA/NGS/$years* ../../../VLBI/DATA/NGS
else
	echo "skipping NGS"
fi
if [ "$vgosDB" = "y" ]; then
	echo "download vgosDB"
	rsync -avzh $user@vievs:/home/members/vievs/VLBI/DATA/vgosDB/$years* ../../../VLBI/DATA/vgosDB
else
	echo "skipping vgosDB"
fi
echo "download VMF1"
scp  $user@vievs:/home/members/vievs/VLBI/TRP/VMF1/y$years.vmf1_r ../../../VLBI/TRP/VMF1
echo "download VMF3"
scp  $user@vievs:/home/members/vievs/VLBI/TRP/VMF3/y$years.vmf3_r ../../../VLBI/TRP/VMF3
echo "download GRAD"
scp  $user@vievs:/home/members/vievs/VLBI/TRP/GRAD/y$years.grad_r ../../../VLBI/TRP/GRAD
echo "download LHG"
scp  $user@vievs:/home/members/vievs/VLBI/TRP/LHG/y$years.lhg_r ../../../VLBI/TRP/LHG
