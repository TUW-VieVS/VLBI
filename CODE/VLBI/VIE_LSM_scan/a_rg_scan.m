
%  07 Jan 2014 by Hana Krasna: estimation of APL regression coefficients

function [A_rg] = a_rg_scan(per_stat,na)

    for j=1:na % all antennas
        for k=1:per_stat(j).total % all observations in one scan
            A_rg(k,j) = per_stat(j).first(k) * per_stat(j).drg(k);    
        end
    end