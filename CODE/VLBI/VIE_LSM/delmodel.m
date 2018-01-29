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
%       'opt'        structure array     (for info. /DOC/opt.doc)
%       'A'          structure array     design matrix of the real-observation equations
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'T'          structure array     estimation intervals for dUT1, xpol, ypol, nutdx, and nutdy
%       'och'        structure array     o-c vector of constraints (zero vector)
%       'n_'         structure array     number of estimates (pwlo or one offset)
%
%   Output:
%       'A'          structure array     design matrix of the real-observation equations
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'T'          structure array     estimation intervals for dUT1, xpol, ypol, nutdx, and nutdy
%       'och'        structure array     o-c vector of constraints (zero vector)
%       'n_'         structure array     number of estimates (pwlo or one offset)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   12 May 2009 by Kamil Teke
%
%   Revision: 
%   19 Aug 2009 by Kamil Teke: options for clock function are added - opt.pw_clk
%   06 Dec 2009 by Kamil Teke: header added
% ************************************************************************ 
function [A,H,Ph,T,och,n_] = delmodel(opt,A,H,Ph,T,och,n_)

if opt.pw_clk == 1  % Deleting rate and quadratic terms from clock function - only pwl offsets (main solution)
    A(2).sm = []; H(2).sm = []; Ph(2).sm = []; och(2).sv = []; 
    for istat = 1 : size(opt.stat,2)
        n_(istat).qclk = 0;
    end
end

if opt.pw_clk == 0  % Deleting pwl offsets + rate + quadratic terms from clock function - (main solution)
    A(1).sm = []; H(1).sm = []; Ph(1).sm = []; och(1).sv = []; 
    A(2).sm = []; H(2).sm = []; Ph(2).sm = []; och(2).sv = []; 
    for istat = 1 : size(opt.stat,2)
        n_(istat).clk = 0;
        n_(istat).qclk = 0;
    end
end

if opt.nutdy.model ~= 1
    A(10).sm = []; H(10).sm = []; Ph(10).sm = []; och(10).sv = []; T.nutdy = [];
 else if opt.dut1.model == 1 && opt.nutdy.constrain == 0
         H(10).sm(1:size(H(10).sm,1),:) = 0; 
         Ph(10).sm(1:size(Ph(10).sm,1),:) = 0;
         och(10).sv = [];
     end
end

if opt.nutdx.model ~= 1
    A(9).sm = []; H(9).sm = []; Ph(9).sm = []; och(9).sv = []; T.nutdx = [];
 else if opt.dut1.model == 1 && opt.nutdx.constrain == 0
         H(9).sm(1:size(H(9).sm,1),:) = 0; 
         Ph(9).sm(1:size(Ph(9).sm,1),:) = 0;
         och(9).sv = [];
     end
end

if opt.dut1.model ~= 1
    A(8).sm = []; H(8).sm = []; Ph(8).sm = []; och(8).sv = []; T.dut1 = [];
 else if opt.dut1.model == 1 && opt.dut1.constrain == 0 
         H(8).sm(1:size(H(8).sm,1),:) = 0; 
         Ph(8).sm(1:size(Ph(8).sm,1),:) = 0;  
         och(8).sv = [];
     end
end

if opt.ypol.model ~= 1 
    A(7).sm = []; H(7).sm = []; Ph(7).sm = []; och(7).sv = []; T.ypol = [];
 else if opt.ypol.model == 1 && opt.ypol.constrain == 0    
         H(7).sm(1:size(H(7).sm,1),:) = 0; 
         Ph(7).sm(1:size(Ph(7).sm,1),:) = 0;
         och(7).sv = [];
     end
end

if opt.xpol.model ~= 1
     A(6).sm = []; H(6).sm = []; Ph(6).sm = []; och(6).sv = []; T.xpol = [];
 else if opt.xpol.model == 1 && opt.xpol.constrain == 0
         H(6).sm(1:size(H(6).sm,1),:) = 0; 
         Ph(6).sm(1:size(Ph(6).sm,1),:) = 0;
         och(6).sv = [];
     end   
end


