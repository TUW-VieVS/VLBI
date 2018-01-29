% ************************************************************************
%
% DESCRIPTION
%   This is a function for get_out_val to get weighted value and its error
%   for repeated EOP at the one time moment.
%	It should exist in the OUT directory. The 
%   function reads the output from get_out_val.m.
%
% CREATED  
%   2015/AUG     M. Madzak, A. Girdiuk
%
% REFERENCES
%
% INPUT
% - values            : two EOP at the one time moment
% - sigmas            : two sigma for EOP correspondingly
%
% OUTPUT
% - weightedVal       : weighted EOP based on sigmas
% - err_of_weightedVal: its error
%
% CHANGES:
%   29-Nov-2016, A.Girdiuk: description is changed wrt introduced new get_out_val function
%                           err_of_weightedVal is returned
% ************************************************************************

function [weightedVal,err_of_weightedVal]=weightedMeanFunction(values,sigmas)

weights     = 1./sigmas.^2;
normWeights = weights/sum(weights);     err_of_weightedVal = sqrt(1./sum(weights));
weightedVal = sum(values.*normWeights);  %   weightedMean

