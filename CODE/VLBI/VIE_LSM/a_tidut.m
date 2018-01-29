% ************************************************************************
%   Description:
%   function to form the design matrix for tidal dUT1 terms
%
%   Reference: 
%
%   Input:										
%      'temp = scan.obs'   structure array   (for info. /DOC/scan.doc)
%      'nobserv'           structure array   total number of observations
%                                            in the session after outliers and bad-quality-coded observations
%                                            eleminated
%
%   Output:
%      'A_tidut'    design matrix (nobserv x 2*ntid)     
%      
%   External calls: 
%   
%   Coded for VieVS: 
%   04 Nov 2011 by Sigrid Boehm
%
%   Revision: 
% ************************************************************************
function [A_tidutc,A_tiduts] = a_tidut(temp,nobserv,c,rad2mas)

ntid = length(temp(1).putc);

A_tidutc = zeros(nobserv,ntid);
A_tiduts = zeros(nobserv,ntid);

for k = 1 : nobserv
    for j = 1 : ntid % loop over ntid columns for tidal terms
        A_tidutc(k,j) = temp(k).putc(j)*c*100*(1/rad2mas);       % Uc amplitudes 
        A_tiduts(k,j) = temp(k).puts(j)*c*100*(1/rad2mas);       % Us amplitudes
    end
end

