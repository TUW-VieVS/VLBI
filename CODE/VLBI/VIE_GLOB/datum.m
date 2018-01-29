% ************************************************************************
%   Description:
%       This function opens the .txt file in '../../DATA/GLOB/TRF/DATUM/'
%       and finds the indices of the stations, which will be used for the
%       NNT/NNR condition.
%       for station coordinates and velocities
%
%
%   Input:
%      file_datum          the chosen file with station names for the
%                          NNT/NNR condition ('../../DATA/GLOB/TRF/DATUM/*.txt')
%      refname             names of antenna, also multiple if there are
%                          breaks in the position, that will be considered
%                          (refnamec/refnamev)
%   Output: 
%      inidant             indices of stations included in the NNT/NNR
%                          condition (wrt refnamec/refnamev)
%      datumant            names of stations included in the NNT/NNR
%                          condition - also multiple
%
%   External calls: 	
%      -               					    											
%       
%   Coded for VieVS: 
%   03 Aug 2010 by Hana Spicakova
%
%   Revision: 
%%


function [inidant,datumant]=datum(file_datum,refname)
    
    inidant=[];
    datumant=[''];
    
    fid=fopen(file_datum);
    if fid~=-1
        datumant=textread(file_datum, '%8c', 'whitespace','','commentstyle','matlab');
        fclose(fid);
    end

    nda=size(datumant,1);
    for i=1:nda
        [inant,xx]=find(strcmp(cellstr(datumant(i,1:8)),cellstr(refname)) == 1);
        inidant=[inidant; inant];
    end
    inidant=unique(inidant); % indices of stations included in the NNT/NNR condition (wrt refname)
    clear datumant
    datumant=refname(inidant,:); % names of stations included in the NNT/NNR condition - also multiple

