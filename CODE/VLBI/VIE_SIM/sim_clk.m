% ************************************************************************
%   Description:
%   This function simulates station clocks as the sum of a random walk and
%   an integrated random walk. The stochastic processes are driven by power
%   spectral densities corresponding to a certain Allan standard deviation
%   (ASD), e.g. 1e-14@50 min, which corresponds to the state-of-the-art
%   clock performance used in VLBI.
% 
%   References: 
%    - Herring et al. 1990 "Geodesy by Radio Interferometry: The
%       Application of Kalman Filtering to the Analysis of Very Long
%       Baseline Interferometry Data"
%    - Böhm et al. 2007 "Simulation of zenith wet delays and clocks"
%
%   Input:										
%      mjd  -  [d] a vector containing the modified julian date of the obs
%      sy1  -  if you want to simulate a clock with an ASD of 1e-14@50 min,
%              sy1 is 1e-14
%      sy2  -  [min] if you want to simulate a clock with an ASD of 1e-14@50 min,
%              sy2 is 50
%      idays   [d] how many days should be simulated
% 
%   Output:
%      a vector of simulated clock values [s] further used in vie_sim
% 
%   External calls: 	
%   ---    
%       
%   Coded for VieVS: 
%   July 2010 by Andrea Pany
%
%   Revision: 
%   Dec 2012 by Sigrid Boehm: enable equal epochs for use with twin telescopes
%
% ************************************************************************

function clk = sim_clk(mjd,sy1,sy2,idays)

% number of observations
num = length(mjd);

% pre-allocation
clk = zeros(num,idays);

% compute time differences between observations in seconds
dts = diff(mjd)*86400;             % [s]

% simulate station clock
cn = 0;
phicr  = sy1^2*sy2*60;                % [sec^2/sec]
phici  = sy1^2/(sy2*60)*3;            % [sec^2/sec^3]
wnsr = sqrt(phicr);                   % [sec/sqrt(sec)]
wnsi = sqrt(phici);                   % [sec/(sec*sqrt(sec))]
ti = zeros(num,1);
tr = zeros(num,1);

for iday = 1:idays
    % random walk
    tr(1) = randn*cn;
    yc = tr(1);
    for i = 2:num
        if dts(i-1)==0
           tr(i) = tr(i-1);
        else
           wnsc = randn*wnsr/sqrt(dts(i-1));        % [sec/sec]
           yc = yc + wnsc*dts(i-1) + randn*cn;      % [sec]
           tr(i) = yc;
        end
    end

    % integrated random walk
    ti(1) = randn*cn;
    yc = ti(1);
    v = 0;
    for i = 2:num
        if dts(i-1)==0
           ti(i) = ti(i-1);
        else
           wnsc = randn*wnsi/sqrt(dts(i-1));        % [sec/sec^2]
           yc = yc + v*dts(i-1) + wnsc/2*dts(i-1)^2 + randn*cn;
           v = v + wnsc*dts(i-1);                   % [sec/sec]
           ti(i) = yc;
        end
    end

    % sum of random walk and integrated random walk
    clk(:,iday) = tr + ti;                               % [sec]
end