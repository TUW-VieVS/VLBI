% ************************************************************************
%   Description:
%   function to form the design matrix of station coordinates global solution: 
%   estimation of antenna axis offset
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
%      'A_ao'    design matrix (nobserv x 1)   design matrix of the station axis offset estimations
%   
%   Coded for VieVS: 
%	04 Oct 2013 by Hana Krasna
%
%   Revision: 
%   10 Oct 2016 by A. Girdiuk: code-optimized: loop removed, 
%							   function assignment is changed in accordance with stwisepar.m
%							   header is added
function [A_ao] = a_AO(obs_per_stat,nobserv,nob)

A_ao(nobserv,1) = 0;

first = obs_per_stat.first; % -1 if it is i1 OR +1 if it is i2 
%total = obs_per_stat.total; % total number of observations that are carried out by the station
%nob = obs_per_stat.nob; % the row numbers of the observations of that station in the oc vector, A, and P matrices 
dAO = obs_per_stat.dAO; % partial derivatives of the observations with respect to AO (one specific station)

% assigning the partial derivatives of station coordinates to the specific rows 

A_ao(nob,1) = first.*dAO;
