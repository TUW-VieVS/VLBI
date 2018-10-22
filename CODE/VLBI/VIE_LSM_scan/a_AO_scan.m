
%  07 Jan 2014 by Hana Krasna: estimation of antenna axis offset

function [A_ao] = a_AO_scan(per_stat,na)

    for j=1:na % all antennas
        for k=1:per_stat(j).total % all observations in one scan
            A_ao(k,j) = per_stat(j).first(k) * per_stat(j).dAO(k);    
        end
    end
