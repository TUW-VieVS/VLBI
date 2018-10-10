% ************************************************************************
%   Description:
%   function to form the design matrix of source coordinates for global solution 
%   (one estimate per session)
%
%   Reference: 
%
%   Input:										
%      'per_source'   structure array   (for info. /DOC/per_source.doc)
%      'nobserv'      matrix (1x1)       total number of observations
%                                        in the session after outliers and bad-quality-coded observations
%                                        eleminated
%
%   Output:
%      'Ara'    design matrix (nobserv x 1)   design matrix of the source's right ascension estimate
%      'Ade'    design matrix (nobserv x 1)   design matrix of the source's declination offset estimate
%      
%   External calls: 	
%   
%   Coded for VieVS: 
%   25 Aug 2009 by Kamil Teke 
%
%   Revision: 
%   07 Oct 2009 by Kamil Teke: header added
%
% ************************************************************************
function [Ara,Ade] = a_source(per_source,nobserv)

Ara(nobserv,1) = 0;
Ade(nobserv,1) = 0;

total = per_source.total; % total number of observations of the specific source in a session
nob = per_source.nob; % the row numbers of the observations to that source in the oc vector, A, and P matrices 
ra = per_source.ra; % [cm/mas] --- partial derivatives of the right ascension coordinates of the source 
de = per_source.de; % [cm/mas] --- partial derivatives of the declination coordinates of the source

 % assigning the partial derivatives to the specific rows of the observations which
 % are to made to the corresponding sourse  
for k = 1 : total
    Ara(nob(k),1) = ra(k);
    Ade(nob(k),1) = de(k); 
end