% ************************************************************************
%   Description:
%   This function simulates white noise values to be added per baseline
%   observation.
%
%   References: 
%   ---
%
%   Input:		
%      wnoise:     standard deviation of the wn process in [ps]
%      nbslobs:    number of baseline observations
%      idays:      how many days should be simulated [d]
% 
%   Output:
%      a vector containing the white noise values in [s]
% 
%   External calls: 
%   ---
%       
%   Coded for VieVS: 
%   July 2010 by Andrea Pany
%
%   Revision: 
%   2016-10-21, A. Hellerschmied: wn preallocated correctly for idays
% ************************************************************************

function wn = sim_wn(wnoise,nbslobs,idays)

wnoise = wnoise*1e-12;      % [s]

wn = zeros(nbslobs,idays);

for iday = 1:idays
    wn(:,iday) = randn(nbslobs,1)*wnoise;
end