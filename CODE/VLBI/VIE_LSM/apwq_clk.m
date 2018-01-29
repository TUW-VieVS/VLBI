% ************************************************************************
%   Description:
%   function to form the design matrix for the clocks as quadratic 
%   functions "Arqclk"  and as piecewise linear ofsets "Apwclk"
%
%   Reference: 
%
%   Input:										
%       'obs_per_stat'        structure array     (for info. /DOC/obs_per_stat.doc)
%       'nobserv'             (1,1)               number of observation in the session
%	    'nob'			vector (1 x n)	   ordered number of observation per station
%       'n_unk'               structure array     number of estimation intervals for clocks, zwd, ngr, egr, xyz    
%       'T_'                  structure array     estimation epochs for clocks, zwd, ngr, egr, xyz
%       'opt'                 structure array     (for info. /DOC/opt.doc)
%
%   Output:
%       'Apwclk'      (nobserv,n_unk.clk+1)          design matrix for pwl clock offsets          
%       'Arqclk'      (noebserv,num. of clocks x 2)  rate and quadratic terms of the clock function
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   12 May 2009 by Kamil Teke
%
%   Revision: 
%   05 Dec 2009 by Kamil Teke: header added
%   10 Oct 2016 by A. Girdiuk: code-optimized: loop removed, 
%							   function assignment is changed in accordance with stwisepar.m
% ************************************************************************ 
function [Apwclk, Arqclk] = apwq_clk(obs_per_stat,nobserv,nob,n_unk,T_,opt)

minute = obs_per_stat.minute;
first = obs_per_stat.first;
no_int_clk = obs_per_stat.no_int_clk;

% number of scans per clock estimate interval [1]

%for inter = 1:n_unk.clk
 %   if inter == n_unk.clk
%        tst = (minute >= T_.clk(inter)).*(minute <= T_.clk(inter+1));
    %else
   %     tst = (minute >= T_.clk(inter)).*(minute < T_.clk(inter+1));
  %  end
 %   stm_clk(inter) = sum(tst); % number of observation to the corresponding source in the 'inter' estimation interval
%end

% Piecewise linear clock offsets
Apwclk(nobserv,n_unk.clk+1) = 0;
% Rate and quadratic terms of the clock function
Arqclk(nobserv,1) = 0;
k = 0; 
for inter = 1:n_unk.clk % number of clock estimate intervals in a session
    for iobs = 1:n_unk.stm_clk(inter) % number of observations of the station in an interval
        k = k + 1;
        
        Apwclk(nob(k),inter) = first(k)*(1-(minute(k)-T_.clk(no_int_clk(k)))/(T_.clk(no_int_clk(k)+1)-T_.clk(no_int_clk(k)))); % [1]
        Apwclk(nob(k),inter+1) = first(k)*(minute(k)-T_.clk(no_int_clk(k)))/(T_.clk(no_int_clk(k)+1)-T_.clk(no_int_clk(k))); % [1]
        
        Arqclk(nob(k),1) = first(k)*(minute(k)-T_.clk(no_int_clk(1)))./60/24; % [days]
        if opt.pw_clk ~= 2
        	Arqclk(nob(k),2) = first(k)*((minute(k)-T_.clk(no_int_clk(1)))./60/24).^2; % [days^2]
        end
    end
end
