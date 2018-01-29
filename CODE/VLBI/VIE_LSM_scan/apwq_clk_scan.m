
% ************************************************************************
%   Description:
%   function to form the design matrix for the clocks as quadratic 
%   functions "Arqclk"  and as piecewise linear ofsets "Apwclk"
%
%   Reference: 
%
%   Input:										
%       'per_stat'   structure array     (for info. /DOC/per_stat.doc)
%       'T_'         structure array      estimation epochs for clocks, zwd, ngr, egr, xyz
%       'opt'        structure array     (for info. /DOC/opt.doc)
%       'na'         (1,1)                Number of antennas
%       't'          (1,num. of scans)    minutes from midnight to the time of the scans
%       'iscan'      (1,1)                Number of the scan that is being analysed
%
%   Output:
%       'Apwclk'     structure array                  design matrix for pwl clock offsets          
%       'Arqclk'     (noebserv,num. of clocks x 2)    rate and quadratic terms of the clock function
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   July 2012 by Claudia Tierno Ros
%
%   Revision: 
%   
% ************************************************************************ 


function [Apwclk, Arqclk,num_inter_qclk] = apwq_clk_scan(per_stat,T_,opt,na,t,iscan)


    for j=1:na
        
        for k=1:per_stat(j).total
              
            Apwclk(j).sm(k,1)= per_stat(j).first(k)*...
                (1-((t(iscan)-T_(j).clk(1))/(T_(j).clk(2)...
                 -T_(j).clk(1)))); 
           
            Apwclk(j).sm(k,2)= per_stat(j).first(k)*...
                ((t(iscan)-T_(j).clk(1))/((T_(j).clk(2))...
                 -T_(j).clk(1))); 
        
            
         
            Arqclk(j).sm(k,1)=per_stat(j).first(k)*(t(iscan)-T_(j).clk_first_obs)/60/24;
      
            if opt.pw_clk~=2
               Arqclk(j).sm(k,2)=per_stat(j).first(k)*((t(iscan)-T_(j).clk_first_obs)/60/24)^2;
            end
        
             
        end
        num_inter_qclk(j)=size(Arqclk(j).sm,2);
    end


