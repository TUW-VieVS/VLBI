% ************************************************************************
%   Description:
%   function to form the design matrix for the shida numbers "a_shida" as 79 estimate
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
%      'A_shida'    design matrix (nobserv x 79)     design matrix of 79 shida numbers for specific frequency bands  
%      
%   External calls: 	
%   
%   Coded for VieVS: 
%   15 Aug 2009 by Kamil Teke and Hana Spicakova
%
%   Revision: 
%   07 Oct 2009 by Kamil Teke: header added
%   23 May 2010 by Hana Spicakova: removed multiplication with -1 of the
%                                  partial derivatives
%   15 Nov 2012 by Hana Krásná: the code was changed to a general way
% ************************************************************************
function [A_shida] = a_shida(temp,nobserv)

nl=length(temp(1).pShida);

A_shida(nobserv,nl) = 0;

for k = 1 : nobserv
    for j = 1 : nl % loop over 79 columns for shida numbers corresponding specific frequency bands
        A_shida(k,j) = temp(k).pShida(j); % assigning observationwise partial  derivatives
    end
end

