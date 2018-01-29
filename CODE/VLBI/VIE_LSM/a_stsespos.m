% ************************************************************************
%   Description:
%   function to form the design matrix of amplitudes of seasonal variations in station positions
%   for global solution 
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
%      	
%   
%   Coded for VieVS: 
%   11 Oct 2011 by Hana Spicakova: based on functions from Kamil Teke
%   10 Oct 2016 by A. Girdiuk: code-optimized: loop removed, 
%							   function assignment is changed in accordance with stwisepar.m
% ************************************************************************
function [AAcr,AAce,AAcn,AAsr,AAse,AAsn] = a_stsespos(obs_per_stat,nobserv,nob)

nt=size(obs_per_stat.pAcr,2);

AAcr(nobserv,nt) = 0;
AAce(nobserv,nt) = 0;
AAcn(nobserv,nt) = 0;
AAsr(nobserv,nt) = 0;
AAse(nobserv,nt) = 0;
AAsn(nobserv,nt) = 0;

first = obs_per_stat.first; % -1 if it is i1 OR +1 if it is i2 
%total = obs_per_stat.total; % total number of observations that are carried out by the station
%nob = obs_per_stat.nob; % the row numbers of the observations of that station in the oc vector, A, and P matrices 


pAcr = obs_per_stat.pAcr; % partial derivatives of the observations with respect to Acr (one specific station)
pAce = obs_per_stat.pAce; % partial derivatives of the observations with respect to Ace (one specific station)
pAcn = obs_per_stat.pAcn; % partial derivatives of the observations with respect to Acn (one specific station)
pAsr = obs_per_stat.pAsr; % partial derivatives of the observations with respect to Asr (one specific station)
pAse = obs_per_stat.pAse; % partial derivatives of the observations with respect to Ase (one specific station)
pAsn = obs_per_stat.pAsn; % partial derivatives of the observations with respect to Asn (one specific station)



% assigning the partial derivatives  to the specific rows 

AAcr(nob,:) = [first(:).*pAcr(:,1) first(:).*pAcr(:,2)];
AAce(nob,:) = [first(:).*pAce(:,1) first(:).*pAce(:,2)];
AAcn(nob,:) = [first(:).*pAcn(:,1) first(:).*pAcn(:,2)];

AAsr(nob,:) = [first(:).*pAsr(:,1) first(:).*pAsr(:,2)];
AAse(nob,:) = [first(:).*pAse(:,1) first(:).*pAse(:,2)];
AAsn(nob,:) = [first(:).*pAsn(:,1) first(:).*pAsn(:,2)];
