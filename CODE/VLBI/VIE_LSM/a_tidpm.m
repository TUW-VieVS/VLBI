% ************************************************************************
%   Description:
%   function to form the design matrix for tidal polar motion terms
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
%      'A_tidpm'    design matrix (nobserv x 2*(ntidp+ntidm))    
%      
%   External calls: 
%   
%   Coded for VieVS: 
%   04 Nov 2011 by Sigrid Boehm
%   
% ************************************************************************
function [A_tidap,A_tidbp,A_tidam,A_tidbm] = a_tidpm(temp,nobserv,c,rad2mas)

ntidp = length(temp(1).ppap);
ntidm = length(temp(1).ppam);

A_tidap = zeros(nobserv,ntidp);
A_tidbp = zeros(nobserv,ntidp);
A_tidam = zeros(nobserv,ntidm);
A_tidbm = zeros(nobserv,ntidm);

for k = 1 : nobserv
    for j = 1 : ntidp % loop over ntidp columns for prograde tidal terms
        A_tidap(k,j) = temp(k).ppap(j)*c*100*(1/rad2mas);       % A+ amplitudes 
        A_tidbp(k,j) = temp(k).ppbp(j)*c*100*(1/rad2mas);       % B+ amplitudes
    end
    
    for l = 1 : ntidm % loop over ntidm columns for retrograde tidal terms
        A_tidam(k,l) = temp(k).ppam(l)*c*100*(1/rad2mas);       % A- amplitudes 
        A_tidbm(k,l) = temp(k).ppbm(l)*c*100*(1/rad2mas);       % B- amplitudes
    end
end

