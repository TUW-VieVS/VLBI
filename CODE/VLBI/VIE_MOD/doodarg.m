% ************************************************************************
%   Description:
%   Computes Doodson's fundamental arguments. Translated from Occam routine
%   matthew.
% 
%   Input:										
%      mjd              Modified Julian Date [d]
%      leap             difference between UTC and TAI in [s]
%                       leap seconds(TAI-UTC)
% 
%   Output:
%      tau              mean lunar time [deg]
%      s                mean tropic longitude of the Moon [deg]
%      h                mean tropic longitude of the Sun [deg]
%      p                mean tropic longitude of the lunar perigee [deg]
%      zns              mean tropic longitude of the ascending lunar node
%                       (decreasing in time)(N') [deg]
%      ps               mean tropic longitude of the perihelion [deg]
% 
%   External calls: 	
%      ---
%       
%   Coded for VieVS: 
%   26 Oct 2008 by Hana Spicakova
%
%   Revision: 
%   19 Jul 2011 by Hana Spicakova: The time (T) was changed to TT
%              (should be TDB, but the difference is negligable)
% ************************************************************************ 

function [tau,s,h,p,zns,ps] = doodarg(mjd,leap)

mjd_ref=51544.5;       % 2000-01-01 at 12UT

T = (mjd - mjd_ref + (leap+32.184)/86400)/36525; 
fhr=(mjd-fix(mjd))*24; % [h]


% equations from
% ftp://tai.bipm.org/iers/conv2010/chapter7/dehanttideinel/STEP2DIU.F
 s   = 218.31664563         + 481267.88194    *T - 0.0014663889*T^2 + ...
        0.00000185139*T^3;
 tau = fhr*15 + 280.4606184 +  36000.7700536  *T + 0.00038793  *T^2 - ...
        0.0000000258 *T^3 - s; % mean lunar time 
 PR  =                             1.396971278*T + 0.000308889 *T^2 + ...
        0.000000021  *T^3 + 0.000000007*T^4;
   s = s+PR; % mean tropic longitude of the Moon
 % mean tropic longitude of the Sun 
 h   = 280.46645    + 36000.7697489    *T + 0.00030322222*T^2 + ...
        0.000000020  *T^3-0.00000000654*T^4;
 % mean tropic longitude of the lunar perigee 
 p   =  83.35324312 +  4069.01363525   *T - 0.01032172222*T^2 - ...
        0.0000124991 *T^3+0.00000005263*T^4;
 % mean tropic longitude of the ascending lunar node (decreasing in
 % time)(N')
 zns = 234.95544499 +  1934.13626197   *T - 0.00207561111*T^2 - ...
        0.00000213944*T^3+0.00000001650*T^4;
 % mean tropic longitude of the perihelion 
 ps  = 282.93734098 +     1.71945766667*T + 0.00045688889*T^2 - ...
        0.00000001778*T^3-0.00000000334*T^4; 
%  tau=reduce_deg(tau);
%  s=reduce_deg(s);
%  h=reduce_deg(h);
%  p=reduce_deg(p);
%  zns=reduce_deg(zns);
%  ps=reduce_deg(ps);
