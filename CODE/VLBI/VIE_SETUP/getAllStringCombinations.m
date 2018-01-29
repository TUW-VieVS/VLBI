% #########################################################################
% #     getAllStringCombinations
% #########################################################################
%
% DESCRITPION
% This function generates all double-combinations from a given set of
%	 strings. Used for baselines from stations, e.g.
%
% AUTHOR 
%   moved in separate function by A. Girdiuk
%
% INPUT
%    {'WETTZELL', 'TIGOCONC', SVETLOE'}
%
% OUTPUT
%  	 {'WETTZELL', 'TIGOCONC';
%	  'WETTZELL', 'SVETLOE';
%	  'TIGOCONC', 'SVETLOE'}
%
% CHANGES
%
function combinations=getAllStringCombinations(inputStrings)

nStat=length(inputStrings);

% preallocate
combinations=cell(sum(1:nStat-1), 2);
rowInd=1;

for k=1:nStat-1
    for m=k+1:nStat
        combinations(rowInd,1)=inputStrings(k);
        combinations(rowInd,2)=inputStrings(m);
        
        rowInd=rowInd+1;
    end
end
