# VieVS-VLBI

VLBI module of VieVS

Written by present and former members of the TU Wien VieVS-team

Administrated by Sigrid Böhm (vievsadmin@tuwien.ac.at)

Documentation: https://vievswiki.geo.tuwien.ac.at

## License
> VieVS-VLBI VLBI module of Vienna VLBI and Satellite Software
>
> Copyright (C) 2019  TU Wien
>
> This program is free software: you can redistribute it and/or modify
> it under the terms of the GNU General Public License as published by
> the Free Software Foundation, either version 3 of the License, or
> (at your option) any later version.
>
> This program is distributed in the hope that it will be useful,
> but WITHOUT ANY WARRANTY; without even the implied warranty of
> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
> GNU General Public License for more details.
>
> You should have received a copy of the GNU General Public License
> along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Reference
If you publish results obtained with VieVS-VLBI, please give credit to the VieVS-team by citing the current reference:

Johannes Böhm, Sigrid Böhm, Janina Boisits, Anastasiia Girdiuk, Jakob Gruber, Andreas Hellerschmied, Hana Krásná, Daniel Landskron, Matthias Madzak, David Mayer, Jamie McCallum, Lucia McCallum, Matthias Schartner, Kamil Teke, Vienna VLBI and Satellite Software (VieVS) for Geodesy and Astrometry, Publications of the Astronomical Society of the Pacific, Vol. 130(986), 044503, 2018. http://iopscience.iop.org/article/10.1088/1538-3873/aaa22b

## Getting started

You need a recent Matlab installation (R2014b or later). VieVS is tested on Windows and Linux.
Download or clone the VLBI repository and(!) the COMMON repository and put them to the same directory.
If your main directory is called 'VieVS' for example, your folder structure should look like this:

    /VieVS/COMMON/
    /VieVS/VLBI/
    
Please be aware that this repository contains more or less only the program code. The data files necessary for processing have to be downloaded from external sources.

### Observation data
* NGS can be downloaded for example from ftp://cddis.nasa.gov/vlbi/ivsdata/ngs/ and have to be stored in year folders (NGS-files have to be uncompressed): `/VLBI/DATA/NGS/yyyy/`
* vgosDB can be downloaded from ftp://cddis.gsfc.nasa.gov/vlbi/ivsdata/vgosdb/ and have to be stored in year folders as you download them (compressed tar files, *.tgz or *.tar.gz): `/VLBI/DATA/vgosDB/yyyy/`
* Master files are needed for the automatic generation of a process list, they can be downloaded from ftp://cddis.gsfc.nasa.gov/pub/vlbi/ivscontrol/ as yearly *.txt files and have to be stored in : `/VLBI/DATA/MASTER`

### Earth Orientation parameters
* EOP C04 can be downloaded from ftp://ftp.iers.org/products/eop/long-term/ . Take for example "eopc04_14_IAU2000.62-now" and store it as `C04_14_1962_now.txt` in `/VLBI/EOP/`.
* EOP finals can be downloaded from ftp://ftp.iers.org/products/eop/rapid/standard/ . Take for example "finals2000A.all" and store it as `finals_all_IAU2000.txt` in `/VLBI/EOP/`.

### Troposphere delay models and ray-tracing data
* Mapping functions: VMF1 and VMF3 can be retrieved from http://vmf.geo.tuwien.ac.at/trop_products/VLBI/ . Choose the yearly files and store them in the respective folder (VMF1 or VMF3) in `/VLBI/TRP/`.
* Gradients are available from http://vmf.geo.tuwien.ac.at/trop_products/VLBI/GRAD/ . Choose the yearly files and store them in the respective folder in `/VLBI/TRP/GRAD`. 
* Ray-traced delays are available from http://vmf.geo.tuwien.ac.at/trop_products/VLBI/RAYTR/RADIATE/. The *.radiate files have to be stored in year folders in `/VLBI/TRP/RAYTRACING_DATA/yyyy/` .

### Running the software
Start Matlab, go to `VLBI/WORK` and type `vievs` in the command window to start the VieVS graphical user interface.
 
