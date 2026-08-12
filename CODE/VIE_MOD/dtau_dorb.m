% ************************************************************************
%   Description:
%   function to compute the partial derivatives of the satellite position
%   w.r.t. the 6 orbital elements and the 9 dynamical parameters
%
%   Reference:
%
%   Input:
%  
%
%   Output:
%       dorb_dt  - partial der. (numerically using time delay tau)
%       dorb_dr  - partial der. (numerically using satellite position)
%       dorb_ana - partial der. (analytically)
%       dorb_frp - partial der. from FRP file for 6 orbit elements (Bernese)
%       dsrp     - partial der. from FRP file for dynamical parameters (Bernese)
%
%   External calls:
%
%   Coded for VieVS:
%   7 May 2025 by H. Wolf 
%
%   Revision:
%
% ************************************************************************
function [dorb_dt, dorb_dr, dorb_ana, dorb_frp, dsrp] = dtau_dorb(GM, scan, sat, dorb, tau, tauC_orb, rad2mas, dsat_gcrf)
    global c       
    tosc = sat.tosc;
    r = scan.crfSat;
    v = scan.v_crfSat;

    %numerical dT
    dorb_dt(1,1) = (tauC_orb(1) - tau) / (dorb(1)) *c ;                % [s/m]*[m/s] -> [] estimates in cm
    dorb_dt(2,1) = (tauC_orb(2) - tau) / (dorb(2)) *c*100;             % [cm]   
    dorb_dt(3,1) = (tauC_orb(3) - tau) / (dorb(3)) *c*100*(1/rad2mas); % [cm/mas]
    dorb_dt(4,1) = (tauC_orb(4) - tau) / (dorb(4)) *c*100*(1/rad2mas); % [cm/mas]
    dorb_dt(5,1) = (tauC_orb(5) - tau) / (dorb(5)) *c*100*(1/rad2mas); % [cm/mas]
    dorb_dt(6,1) = (tauC_orb(6) - tau) / (dorb(6)) *c*100*(1/rad2mas); % [cm/mas]

    %analytically
    [a,e,i,Omega,argp,t00] = xyzele(GM,scan.mjd,r,v);
    [drdorb_ana] = rpartn(GM, scan.mjd, sat.tosc, a, e, i, Omega, argp, t00);

    %numerical dR
    [drdorb_dr] = drdorb_num(GM, scan.mjd, tosc, r, v);  %numerical drdpar

    %FRP
    dorb_frp = zeros(6,3);
    dsrp = zeros(9,3);

    %transform units
    dorb_dr(1,1) = dot(dsat_gcrf',drdorb_dr(1,:)');                    % [m/m]     -> [] estimates in cm
    dorb_dr(2,1) = dot(dsat_gcrf',drdorb_dr(2,:)') *100;               % [cm]      -> estimates dimensionless
    dorb_dr(3,1) = dot(dsat_gcrf',drdorb_dr(3,:)') *100*(1/rad2mas);   % [cm/mas]  -> estimates in mas
    dorb_dr(4,1) = dot(dsat_gcrf',drdorb_dr(4,:)') *100*(1/rad2mas);   % [cm/mas]  -> estimates in mas
    dorb_dr(5,1) = dot(dsat_gcrf',drdorb_dr(5,:)') *100*(1/rad2mas) ;  % [cm/mas]  -> estimates in mas
    dorb_dr(6,1) = dot(dsat_gcrf',drdorb_dr(6,:)') *100*(1/rad2mas);   % [cm/mas]  -> estimates in mas 

    dorb_ana(1,1) = dot(dsat_gcrf',drdorb_ana(1,:)');                    % [m/m]     -> [] estimates in cm
    dorb_ana(2,1) = dot(dsat_gcrf',drdorb_ana(2,:)') *100;               % [cm]      -> estimates dimensionless
    dorb_ana(3,1) = dot(dsat_gcrf',drdorb_ana(3,:)') *100*(1/rad2mas);   % [cm/mas]  -> estimates in mas
    dorb_ana(4,1) = dot(dsat_gcrf',drdorb_ana(4,:)') *100*(1/rad2mas);   % [cm/mas]  -> estimates in mas
    dorb_ana(5,1) = dot(dsat_gcrf',drdorb_ana(5,:)') *100*(1/rad2mas) ;  % [cm/mas]  -> estimates in mas
    dorb_ana(6,1) = dot(dsat_gcrf',drdorb_ana(6,:)') *100*(1/rad2mas);   % [cm/mas]  -> estimates in mas 

end