% ************************************************************************
%   Description:
%   function to form the design matrix of source coordinates for global solution 
%   (one estimate per session)
%
%   Reference: 
%
%   Input:										
%      'per_source'   structure array   (for info. /DOC/per_source.doc)
%
%   Output:
%      'Ara'    design matrix (nobserv per scan x 1)   design matrix of the source's right ascension estimate
%      'Ade'    design matrix (nobserv per scan x 1)   design matrix of the source's declination offset estimate
%      
%   External calls: 	
%   
%   Coded for VieVS: 
%   July 2012 by Claudia Tierno Ros 
%
%   Revision: 
%
% ************************************************************************

function [Ara,Ade] = a_source_scan(per_source)

  
      for k=1:per_source.total
            Ara(k,1)=per_source.ra(k);
            Ade(k,1)=per_source.de(k);
      end
    

