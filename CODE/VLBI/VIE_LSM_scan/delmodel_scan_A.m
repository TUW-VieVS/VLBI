
% ************************************************************************
%   Description:
%   function to exclude models, or exclude soft constraints 
%   excluded models OR soft constraints excluded models: 
%   pwl offsets + rate + quadratic terms from clock function - (main solution)
%   nutation dx (~dEPS), nutation dy (~dPSI)  ---  CELESTIAL POLE OFFSETS
%   XPOL, YPOL, dUT1 --- POLAR MOTION 
%
%   Reference: 
%
%   Input:	
%       'opt'             structure array       (for info. /DOC/opt.doc)
%       'A_session'       structure array        design matrix of the real-observation equations
%
%   Output:
%       'A_session'     structure array     design matrix of the real-observation equations
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   31 July 2012 by Claudia Tierno Ros
%
%   Revision: 
%
% ************************************************************************

function [A_session] = delmodel_scan_A(opt,A_session)

if opt.pw_clk == 1  % Deleting rate and quadratic terms from clock function - only pwl offsets (main solution)
    A_session(2).sm = [];
    
end

if opt.pw_clk == 0  % Deleting pwl offsets + rate + quadratic terms from clock function - (main solution)
    A_session(1).sm = [];  
    A_session(2).sm = []; 
    
end

if opt.nutdy.model ~= 1
    A_session(10).sm = []; 
 
end

if opt.nutdx.model ~= 1
    A_session(9).sm = []; 
 
end

if opt.dut1.model ~= 1
    A_session(8).sm = []; 
 
end

if opt.ypol.model ~= 1 
    A_session(7).sm = [];
 
end

if opt.xpol.model ~= 1
     A_session(6).sm = []; 
   
end