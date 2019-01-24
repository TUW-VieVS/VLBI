# VieVS-VLBI

VLBI module of VieVS

Written by present and former members of the TU Wien VieVS-team

Administrated by Sigrid Böhm (sigrid.boehm@geo.tuwien.ac.at)

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
    
Please be aware that this repository containes more or less only the program code. The data files necessary for processing have to be downloaded from external sources.

### Observation data
* NGS can be downloaded for example from ftp://cddis.nasa.gov/vlbi/ivsdata/ngs/ and have to be stored in year folders (NGS-files have to be uncompressed): `/VLBI/DATA/NGS/yyyy/`
* vgosDB can be downloaded from https://ivscc.gsfc.nasa.gov/hfeop_wg/vgosDB_IVS/ and have to be stored in year folders as you download them (compressed tar files, *.tgz or *.tar.gz): `/VLBI/DATA/vgosDB/yyyy/`
 
