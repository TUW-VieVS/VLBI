% ************************************************************************
%   Description:
%   function to form the design matrix for the source coordinates as
%   piecewise linear ofset functions
%
%   Reference: 
%
%   Input:										
%       'per_source'        structure array     (for info. /DOC/per_source.doc)
%       'T_source'          (1,2)               Estimation epochs for sources coord.
%
%   Output:
%       'Apw_ra'      (nobserv/scan,number of right ascension pwlo unknowns)  design matrix for ra source coord.           
%       'Apw_de'      (nobserv/scan,number of declination pwlo unknowns)      design matrix for de source coord. 
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   July 2012 by Claudia Tierno Ros
%
%   Revision: 
% 
% ************************************************************************ 


function [Apw_ra,Apw_de] = apw_source_scan(per_source,T_source)

 
      for k=1:per_source.total
           
           Apw_ra(k,1) = (1-(per_source.minute-T_source(1))/...
            (T_source(2)-T_source(1)))*per_source.ra(k);

           Apw_ra(k,2) = ((per_source.minute-T_source(1))/...
            (T_source(2)-T_source(1)))*per_source.ra(k);
        
           Apw_de(k,1) = (1-(per_source.minute-T_source(1))/...
            (T_source(2)-T_source(1)))*per_source.de(k);

           Apw_de(k,2) = ((per_source.minute-T_source(1))/...
            (T_source(2)-T_source(1)))*per_source.de(k);
       end
    
