% ************************************************************************
%   Description:
%   function to form the design matrix of station coordinates for single and global solution 
%   (one estimate per session)
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
%      'Ax'    design matrix (nobserv x 1)   design matrix of the station dx coordinate offset estimations
%      'Ay'    design matrix (nobserv x 1)   design matrix of the station dy coordinate offset estimations
%      'Az'    design matrix (nobserv x 1)   design matrix of the station dz coordinate offset estimations
%      
%   External calls: 	
%   
%   Coded for VieVS: 
%   15 Aug 2009 by Kamil Teke 
%
%   Revision: 
%   07 Oct 2009 by Kamil Teke: header added
%   10 Oct 2016 by A. Girdiuk: code-optimized: loop removed, 
%							   function assignment is changed in accordance with stwisepar.m
%
% ************************************************************************
function [Ax,Ay,Az] = a_xyz(obs_per_stat,nobserv,nob)

Ax(nobserv,1) = 0;
Ay(nobserv,1) = 0;
Az(nobserv,1) = 0;

first = obs_per_stat.first; % -1 if it is i1 OR +1 if it is i2 
%total = obs_per_stat.total; % total number of observations that are carried out by the station
%nob = obs_per_stat.nob; % the row numbers of the observations of that station in the oc vector, A, and P matrices 
dx = obs_per_stat.dx; % partial derivatives of the observations with respect to dx (one specific station)
dy = obs_per_stat.dy; % partial derivatives of the observations with respect to dy (one specific station)
dz = obs_per_stat.dz; % partial derivatives of the observations with respect to dz (one specific station)

% assigning the partial derivatives of station coordinates to the specific rows 

Ax(nob,1) = first.*dx;
Ay(nob,1) = first.*dy;
Az(nob,1) = first.*dz;  
