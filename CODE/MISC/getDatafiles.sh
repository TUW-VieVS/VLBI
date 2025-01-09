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

echo "download NTSL from IMLS? (y for yes)"
read imls

echo "download EOPs"
scp $user@hg-proc:/data/DATA4VLBI/EOP/*.txt ../../EOP/

echo "download NTAL"
scp $user@hg-proc:/data/DATA4VLBI/NTSL/NTAL/VIE/vie_y$years.ntal_r ../../NTSL/NTAL/VIE

if [ "$imls" = "y" ]; then
	echo "download NTSL/..../IMLS"
	scp $user@hg-proc:/data/DATA4VLBI/NTSL/NTAL/IMLS/imls_y$years.ntal_r ../../NTSL/NTAL/IMLS
    scp $user@hg-proc:/data/DATA4VLBI/NTSL/NTOL/IMLS/imls_y$years.ntol_r ../../NTSL/NTOL/IMLS
    scp $user@hg-proc:/data/DATA4VLBI/NTSL/HYDL/IMLS/imls_y$years.hydl_r ../../NTSL/HYDL/IMLS
else
	echo "skipping NTSL/..../IMLS"
fi

if [ "$ngs" = "y" ]; then
	echo "download NGS"
	rsync -avzh $user@hg-proc:/data/DATA4VLBI/DATA/NGS/$years* ../../DATA/NGS
else
	echo "skipping NGS"
fi
if [ "$vgosDB" = "y" ]; then
	echo "download vgosDB"
	rsync -avzh $user@hg-proc:/data/DATA4VLBI/DATA/vgosDB/$years* ../../DATA/vgosDB
else
	echo "skipping vgosDB"
fi
if [ "$radiate" = "y" ]; then
        echo "download RAYTRACING_DATA"
        rsync -avzh $user@hg-proc:/data/DATA4VLBI/TRP/RAYTRACING_DATA/$years* ../../TRP/RAYTRACING_DATA
else
        echo "skipping RAYTRACING_DATA"
fi
echo "download VMF1"
scp  $user@hg-proc:/data/DATA4VLBI/TRP/VMF1/y$years.vmf1_r ../../TRP/VMF1
echo "download VMF3"
scp  $user@hg-proc:/data/DATA4VLBI/TRP/VMF3/y$years.vmf3_r ../../TRP/VMF3
echo "download GRAD"
scp  $user@hg-proc:/data/DATA4VLBI/TRP/GRAD/y$years.grad_r ../../TRP/GRAD
