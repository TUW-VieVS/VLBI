% ************************************************************************
%   Description:
%   function that creates t and T so that they can be used in function splitx 
%
%   Reference: 
%
%   Input:										
%       'opt'            structure array            (for info. /DOC/opt.doc)
%       't_last'         (1,num. of scans)          minutes from midnight to the last scan of the antenna
%       't_first'        (1,num. of antennas)       minutes from midnight to the 1st scan of the antenna
%       'na'             (1,1)                      number of antennas
%       'ns'             (1,1)                      number of sources
%       'num_inter_XXX'  (1,num. antennas)/(1,1)    number of interpolation intervals of XXX
%       'nistat'         (1,1)                      number of the reference station
%
%   Output:
%       't'                 structure array         estimation intervals of clk, zwd, ngr, egr, xyz
%       'T'                 structure array         estimation intervals for dUT1, xpol, ypol, nutdx, and nutdy
%       'tso'               structure array         estimation intervals of source coor. for each source 
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   29 August 2012 by Claudia Tierno Ros
%
%   Revision: 
%   
% ************************************************************************




function [t T tso]=tTscan2lsm(opt, t_first,t_last,na,nistat,num_inter_xpol,num_inter_ypol,num_inter_dut1,num_inter_nutdx,num_inter_nutdy,num_inter_sou,ns)


for i=1:na
    
    % Estimation intervals for clock
    if i~=nistat
       t1 = floor(t_first(i)/opt.int_clk)*opt.int_clk;
       t2 = ceil(t_last(i)/opt.int_clk)*opt.int_clk;
       t(i).clk=t1:opt.int_clk:t2;
    else t(i).clk=0;
    end
    
    % Estimation intervals for zwd
    if opt.stat(i).zwd_inc==1
     t1 = floor(t_first(i)/opt.stat(i).int_zwd)*opt.stat(i).int_zwd;
     t2 = ceil(t_last(i)/opt.stat(i).int_zwd)*opt.stat(i).int_zwd;
     t(i).zwd = t1:opt.stat(i).int_zwd:t2; 
    else
        t(i).zwd=0;
    end
    
    
     % Estimation intervals for ngr
     if opt.stat(i).ngr_inc==1
     t1 = floor(t_first(i)/opt.stat(i).int_ngr)*opt.stat(i).int_ngr;
     t2 = ceil(t_last(i)/opt.stat(i).int_ngr)*opt.stat(i).int_ngr;
     t(i).ngr = t1:opt.stat(i).int_ngr:t2; 
     else
         t(i).ngr=0;
     end
     
 
     % Estimation intervals for egr
     if opt.stat(i).egr_inc==1
     t1 = floor(t_first(i)/opt.stat(i).int_egr)*opt.stat(i).int_egr;
     t2 = ceil(t_last(i)/opt.stat(i).int_egr)*opt.stat(i).int_egr;
     t(i).egr = t1:opt.stat(i).int_egr:t2; 
     else
         t(i).egr=0;
     end
     
     % Estimation intervals for xyz
     if opt.stat(i).xyz_inc==1
     t1 = floor(t_first(i)/opt.stat(i).int_xyz)*opt.stat(i).int_xyz;
     t2 = ceil(t_last(i)/opt.stat(i).int_xyz)*opt.stat(i).int_xyz;
     t(i).xyz = t1:opt.stat(i).int_xyz:t2; 
     else
         t(i).xyz=0;
     end
    
     
     
    
end

if opt.xpol.model==1
T.xpol=0:opt.xpol.int:num_inter_xpol*opt.xpol.int;
else
    T.xpol=0;
end

if opt.ypol.model==1
T.ypol=0:opt.ypol.int:num_inter_ypol*opt.ypol.int;
else
    T.ypol=0;
end

if opt.dut1.model==1
T.dut1=0:opt.dut1.int:num_inter_dut1*opt.dut1.int;
else
    T.dut1=0;
end

if opt.nutdx.model==1
T.nutdx=0:opt.nutdx.int:num_inter_nutdx*opt.nutdx.int;
else
    T.nutdx=0;
end

if opt.nutdy.model==1
T.nutdy=0:opt.nutdy.int:num_inter_nutdy*opt.nutdy.int;
else
    T.nutdy=0;
end



for i=1:ns
    tso(i).sources=0:opt.sour_int_rade:num_inter_sou(i)*opt.sour_int_rade;
    
end





