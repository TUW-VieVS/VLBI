% ************************************************************************
%   Description:
%   Finds the correct external ionospheric delay from the iondata file.
% 
%   References: 
%
%   Input:										
% 
%   Output:
% 
%   External calls: 	   
%       
%   Coded for VieVS: 
%   22 May 2012 by Lucia Plank
%
%   Revision: 
%    2016-10-17, A. Hellerschmied: Code revised for VieVS 3.0
%   
% ************************************************************************
function [dion1,dion2] = get_iondel(iondata, tim, statNameI1, statNameI2)
    % get lines with same date as current scan mjd
    ionLineIndDate = iondata{3}==tim(1) & iondata{4}==tim(2) & ...
                     iondata{5}==tim(3) & iondata{6}==tim(4) & ...
                     iondata{7}==tim(5) & iondata{8}==tim(6);
                 

    % find lines also containing the correct station
    ionLineInd1 = strcmp(iondata{9}, deblank(statNameI1)) & ionLineIndDate;
    ionLineInd2 = strcmp(iondata{9}, deblank(statNameI2)) & ionLineIndDate;

    % see if there is one line for iono delays of both station
    dion1 = 0; 
    dion2 = 0;
    if sum(ionLineInd1) ~= 1 
        fprintf('No (or more than one) trp line found for station %s and time %4d-%02d-%02d %02d:%02d:%05.2f  \n\n', statNameI1, tim(1), tim(2), tim(3), tim(4), tim(5), tim(6));
    elseif sum(ionLineInd2)~=1
        fprintf('No (or more than one) trp line found for station %s and time %4d-%02d-%02d %02d:%02d:%05.2f  \n\n', statNameI2, tim(1), tim(2), tim(3), tim(4), tim(5), tim(6));
    else % there is one line found for both iono delays
        % 2-1... apply to obs
        dion1 = iondata{14}(ionLineInd1);
        dion2 = iondata{14}(ionLineInd2);
    end
        
        