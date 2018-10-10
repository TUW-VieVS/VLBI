% ************************************************************************
%   Description:
%   function to form the design matrix for the love numbers "a_love" as 77
%   estimate
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
%      'A_love'    design matrix (nobserv x 77)     design matrix of 77 love numbers for specific frequency bands  
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
function [A_love] = a_love(temp,nobserv)

nl=length(temp(1).pLove);

A_love(nobserv,nl) = 0;

for k = 1 : nobserv
    for j = 1 : nl % loop over 77 columns for love numbers corresponding specific frequency bands
        A_love(k,j) = temp(k).pLove(j); % assigning observationwise partial derivatives
    end
end