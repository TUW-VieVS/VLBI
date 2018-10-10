% ************************************************************************
%   Description:
%   function to form the design matrix of station coordinates global solution: 
%    estimation of APL regression coefficients
%
%   Reference: 
%
%   Input:										
%      'obs_per_stat'   structure array   (for info. /DOC/obs_per_stat.doc)
%      'nobserv'        matrix (1x1)       total number of observations
%                                          in the session after outliers and bad-quality-coded observations
%                                          eleminated
%	   'nob'			vector (1 x n)	   ordered number of observation per station
%
%   Output:
%      'A_ao'    design matrix (nobserv x 1)   design matrix of APL regression coefficients
%   
%   Coded for VieVS: 
%	05 Dec 2013 by Hana Krasna
%
%   Revision: 
%   10 Oct 2016 by A. Girdiuk: code-optimized: loop removed, 
%							   function assignment is changed in accordance with stwisepar.m
%							   header is added

function [AA_rg] = a_rg(obs_per_stat,nobserv,nob)

AA_rg(nobserv,1) = 0;

first = obs_per_stat.first; % -1 if it is i1 OR +1 if it is i2 
%total = obs_per_stat.total; % total number of observations that are carried out by the station
%nob = obs_per_stat.nob; % the row numbers of the observations of that station in the oc vector, A, and P matrices 
drg = obs_per_stat.drg; % partial derivatives of the observations with respect to RG (one specific station)

% assigning the partial derivatives of station coordinates to the specific rows 

AA_rg(nob,1) = first.*drg;
