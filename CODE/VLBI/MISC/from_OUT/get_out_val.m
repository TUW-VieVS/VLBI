% ************************************************************************
%
% DESCRIPTION
%   This is a function for eop_out to get weighted value 
%   for repeated EOP at the one time moment.
%	It should exist in the OUT directory. The 
%   function reads the output from eop_out.m.
%
% CREATED  
%   2016/Nov     A. Girdiuk
%
% REFERENCES
%
% INPUT
% - values            : two EOP at the one time moment
% - sigmas            : two sigma for EOP correspondingly
%
% OUTPUT
% - weightedVal       : weighted EOP based on sigmas
%
%   External calls:
%       /OUT/weightedMeanFunction.m
%
% CHANGES:
%
% ************************************************************************

function [out_val,out_err] = get_out_val(val_keep,val_e_keep)

if sum(isnan(val_keep))==0
    [out_val,out_err]  =weightedMeanFunction(val_keep  ,val_e_keep);
else
    index = ~isnan(val_keep);
    if sum(index)==1
        out_val  =val_keep(index);  out_err=val_e_keep(index);
    elseif sum(~index)~=length(index)
        [out_val,out_err]  =weightedMeanFunction(val_keep(index),val_e_keep  (index));
    else
        out_val = nan;  out_err = nan;
    end
end
