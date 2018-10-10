% ************************************************************************
%   Description:
%   function to form the design matrix for the source coordinates as
%   piecewise linear ofset functions
%
%   Reference: 
%
%   Input:										
%       'per_source'        structure array     (for info. /DOC/per_source.doc)
%       'nobserv'           (1,1)               number of observation in the session
%       'n_unk'             structure array     number of estimation intervals for clocks, zwd, ngr, 
%                                               egr, stations coor. , sources coord.    
%       'T_'                structure array     estimation epochs for clocks, zwd, ngr, 
%                                               egr, stations coor. , sources coord.
%
%   Output:
%       'Apw_ra'      (nobserv,number of right ascension pwlo unknowns)  design matrix for ra source coord.           
%       'Apw_de'      (nobserv,number of declination pwlo unknowns)      design matrix for de source coord. 
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   25 Aug 2009 by Kamil Teke
%
%   Revision: 
%   06 Dec 2009 by Kamil Teke: header added
% ************************************************************************ 
function [Apw_ra,Apw_de] = apw_source(per_source,nobserv,n_unk,T_)

Apw_ra(nobserv,n_unk.sou+1) = 0;
Apw_de(nobserv,n_unk.sou+1) = 0;

minute = per_source.minute;
no_int_rade = per_source.no_int_rade;
total = per_source.total;
nob = per_source.nob;
ra = per_source.ra; % [cm/mas]
de = per_source.de; % [cm/mas]

% number of scans per source coor. estimate interval [1]
m = []; 
for inter = 1:n_unk.sou 
    n = 0; 
    for itim = 1:total 
        if inter == n_unk.sou
            tst= (minute(itim) >= T_.source(inter) && minute(itim) <= T_.source(inter+1));
        else
            tst= (minute(itim) >= T_.source(inter) && minute(itim) < T_.source(inter+1));
        end
        if tst
           n = n + 1; 
           stm_sou(inter) = n; % number of observation to the corresponding source in the 'inter' estimation interval
        else 
        end  
    end 
end

k = 0;
for inter = 1:n_unk.sou % number of source coordinate estimate intervals in a session
    for iobs = 1:stm_sou(inter) % number of obs. of the source in an estimation interval
        k = k + 1;
        Apw_ra(nob(k),inter) = (1-(minute(k)-T_.source(no_int_rade(k)))/...
            (T_.source(no_int_rade(k)+1)-T_.source(no_int_rade(k))))*ra(k);
        Apw_ra(nob(k),inter+1) = ((minute(k)-T_.source(no_int_rade(k)))/...
            (T_.source(no_int_rade(k)+1)-T_.source(no_int_rade(k))))*ra(k);
        
        Apw_de(nob(k),inter) = (1-(minute(k)-T_.source(no_int_rade(k)))/...
            (T_.source(no_int_rade(k)+1)-T_.source(no_int_rade(k))))*de(k);
        Apw_de(nob(k),inter+1) = ((minute(k)-T_.source(no_int_rade(k)))/...
            (T_.source(no_int_rade(k)+1)-T_.source(no_int_rade(k))))*de(k);
    end
end