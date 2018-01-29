% ************************************************************************
%   Description:
%   function that creates n_ so that it can be used in function splitx 
%
%   Reference: 
%
%   Input:										
%       'na'             (1,1)                      number of antennas
%       'num_inter_XXX'  (1,num. antennas)/(1,1)    number of interpolation intervals of XXX
%       'nistat'         (1,1)                      number of the reference station
%
%   Output:
%       'n_'                structure array         number of estimates (pwlo or one offset)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   29 August 2012 by Claudia Tierno Ros 
%
%   Revision: 
%   
% ************************************************************************


function [n_]=nscan2lsm(num_inter_clk,num_inter_egr,num_inter_ngr,num_inter_zwd,na,nistat,num_inter_qclk,opt)


for i=1:na
   
    if i~=nistat
        n_(i).clk=num_inter_clk(i)+1;
        n_(i).qclk=num_inter_qclk(i);
    else n_(i).clk=0;
        n_(i).qclk=0;
    end
   
    if opt.stat(i).zwd_inc==1
    n_(i).zwd=num_inter_zwd(i)+1;
    else
        n_(i).zwd=0;
    end
    
    if opt.stat(i).ngr_inc==1
    n_(i).ngr=num_inter_ngr(i)+1;
    else
         n_(i).ngr=0;
    end
        
    if opt.stat(i).egr_inc==1
    n_(i).egr=num_inter_egr(i)+1;
    else
    n_(i).egr=0;
    end
    
    if opt.stat(i).xyz_inc==1
    n_(i).xyz=1;
    else
    n_(i).xyz=0;
    end
    
    
end





