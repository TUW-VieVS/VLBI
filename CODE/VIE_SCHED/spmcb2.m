% #########################################################################
% #     spmcb2.m
% #########################################################################
%
% DESCRITPION
%   Calculate the possible grouping. 
%
% AUTHOR 
%   SUN Jing (2010-04-06)
%
% INPUT
%   downstanum1     
%   downstasn1      
%   downstanum2    
%   downstasn2      
%
% OUTPUT
%   np
%   pmcb
%
% COUPLING
%
%
% CHANGES
%   2015-03-02, A. Hellerschmied: Use function "nchoosek" instead of "combntns", 
%       because the latter one won't be available in the next MATLAB release.
%       Header added.

function [np, pmcb] = spmcb2(downstanum1, downstasn1, downstanum2, downstasn2)

% check the common station
staln = 0;
for ista = 1 : downstanum1
    staid1 = downstasn1(ista); 
    if (size(find(downstasn2(1:downstanum2)==staid1),2) == 1)
        staln = staln + 1;
        stals(staln) = staid1;
    end
end
downstanum11   = 0;
downstasn11(1) =0;
for ista = 1 : downstanum1
    staid1 = downstasn1(ista); 
    if (size(find(downstasn2(1:downstanum2)==staid1),2) == 0)
        downstanum11 = downstanum11 + 1;
        downstasn11(downstanum11) = staid1;
    end
end
downstanum22   = 0;
downstasn22(1) =0;
for ista = 1 : downstanum2
    staid2 = downstasn2(ista); 
    if (size(find(downstasn1(1:downstanum1)==staid2),2) == 0)
        downstanum22 = downstanum22 + 1;
        downstasn22(downstanum22) = staid2;
    end
end

% possibility number
np = 0;

equal = 0;
not_equal = 0;
    
for il = 1 : staln
    % tmpm = combntns(stals, il); % Function "combntns" not available any more after MATLAB2014b
    tmpm = nchoosek(stals, il);
    
    for im = 1 : size(tmpm,1)
        np = np + 1;
        pmcb(np).downstanum1 = downstanum11 + il;
        pmcb(np).downstasn1(1:downstanum11) = downstasn11(1:downstanum11);
        pmcb(np).downstasn1(downstanum11+1:downstanum11+il) = tmpm(im,1:il);
        sl = 0;
        for i = 1 : staln
            if (size(find(tmpm(im,1:il)==stals(i)),2) == 0)
                sl = sl + 1;
                ss(sl) = stals(i);
            end
        end
        pmcb(np).downstanum2 = downstanum22 + sl;
        pmcb(np).downstasn2(1:downstanum22)  = downstasn22(1:downstanum22);
        for is = 1 : sl
            pmcb(np).downstasn2(downstanum22+is) = ss(is);
        end
    end
end
np = np + 1;
pmcb(np).downstanum1 = downstanum11;
pmcb(np).downstasn1(1:downstanum11) = downstasn11(1:downstanum11);
pmcb(np).downstanum2 = downstanum22 + staln;
pmcb(np).downstasn2(1:downstanum22) = downstasn22(1:downstanum22);
if (staln > 0)
    pmcb(np).downstasn2(downstanum22+1:downstanum22+staln) = stals(1:staln);
end


