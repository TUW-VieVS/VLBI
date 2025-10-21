% ************************************************************************
%   Description:
%   function to compute the partial derivatives of the satellite position
%   w.r.t. the 6 keplerian elements.
%
%   Reference:
%
%   Input:
%  
%
%   Output:
%       pdKepEle_dT  - partial der. (numerically using time delay tau)
%       pdKepEle_dR  - partial der. (numerically using satellite position)
%       pdKepEle_ana - partial der. (analytically)
%       pdKepEle_FRP - partial der. from FRP file (Bernese)
%
%   External calls:
%
%   Coded for VieVS:
%   7 May 2025 by H. Wolf 
%
%   Revision:
%
% ************************************************************************
function [pdKepEle_dT, pdKepEle_dR, pdKepEle_ana, pdKepEle_FRP] = compute_pd_KepEle(GM, scan, sat, satnum, dKepEle, tau, tauC_KepEle1, tauC_KepEle2, tauC_KepEle3, tauC_KepEle4, tauC_KepEle5, tauC_KepEle6, rad2mas, coerpr,trpr,hrpr,tbnd,sclpar, pdSatPosGCRF, parameter)
    global c       
    tosc = sat.tosc;
    r = scan.crfSat;
    v = scan.v_crfSat;

    %numerical dT
    pdKepEle_dT(1,1) = (tauC_KepEle1 - tau) / (dKepEle(1)) *c ;                % [s/m]*[m/s] -> [] estimates in cm
    pdKepEle_dT(2,1) = (tauC_KepEle2 - tau) / (dKepEle(2)) *c*100;             % [cm]   
    pdKepEle_dT(3,1) = (tauC_KepEle3 - tau) / (dKepEle(3)) *c*100*(1/rad2mas); % [cm/mas]
    pdKepEle_dT(4,1) = (tauC_KepEle4 - tau) / (dKepEle(4)) *c*100*(1/rad2mas); % [cm/mas]
    pdKepEle_dT(5,1) = (tauC_KepEle5 - tau) / (dKepEle(5)) *c*100*(1/rad2mas); % [cm/mas]
    pdKepEle_dT(6,1) = (tauC_KepEle6 - tau) / (dKepEle(6)) *c*100*(1/rad2mas); % [cm/mas]

    %analytically
    [a,e,i,Omega,argp,t00] = xyzele(GM,scan.mjd,r,v);
    [drdpar_ana] = rpartn(GM, scan.mjd, sat.tosc, a, e, i, Omega, argp, t00);

    %numerical dR
    [drdpar_dR] = drdorb_num(GM, scan.mjd, tosc, r, v);  %numerical drdpar

    %FRP
    if parameter.lsmopt.KepEle.estKepEle_FRP && parameter.lsmopt.KepEle.estKepEle
        [drdpar_FRP,~] = getrpr(sat.fso_name,scan.time_gps,coerpr,trpr,hrpr,tbnd,sclpar,satnum);
    else
        drdpar_FRP = zeros(6,3);
        %if parameter.lsmopt.KepEle.estKepEle_FRP == 1
        %    fprintf('    - Error: No FRP File!')
        %end
    end

    %transform units
    pdKepEle_dR(1,1) = dot(pdSatPosGCRF',drdpar_dR(1,:)');                    % [m/m]     -> [] estimates in cm
    pdKepEle_dR(2,1) = dot(pdSatPosGCRF',drdpar_dR(2,:)') *100;               % [cm]      -> estimates dimensionless
    pdKepEle_dR(3,1) = dot(pdSatPosGCRF',drdpar_dR(3,:)') *100*(1/rad2mas);   % [cm/mas]  -> estimates in mas
    pdKepEle_dR(4,1) = dot(pdSatPosGCRF',drdpar_dR(4,:)') *100*(1/rad2mas);   % [cm/mas]  -> estimates in mas
    pdKepEle_dR(5,1) = dot(pdSatPosGCRF',drdpar_dR(5,:)') *100*(1/rad2mas) ;  % [cm/mas]  -> estimates in mas
    pdKepEle_dR(6,1) = dot(pdSatPosGCRF',drdpar_dR(6,:)') *100*(1/rad2mas);   % [cm/mas]  -> estimates in mas 

    pdKepEle_ana(1,1) = dot(pdSatPosGCRF',drdpar_ana(1,:)');                    % [m/m]     -> [] estimates in cm
    pdKepEle_ana(2,1) = dot(pdSatPosGCRF',drdpar_ana(2,:)') *100;               % [cm]      -> estimates dimensionless
    pdKepEle_ana(3,1) = dot(pdSatPosGCRF',drdpar_ana(3,:)') *100*(1/rad2mas);   % [cm/mas]  -> estimates in mas
    pdKepEle_ana(4,1) = dot(pdSatPosGCRF',drdpar_ana(4,:)') *100*(1/rad2mas);   % [cm/mas]  -> estimates in mas
    pdKepEle_ana(5,1) = dot(pdSatPosGCRF',drdpar_ana(5,:)') *100*(1/rad2mas) ;  % [cm/mas]  -> estimates in mas
    pdKepEle_ana(6,1) = dot(pdSatPosGCRF',drdpar_ana(6,:)') *100*(1/rad2mas);   % [cm/mas]  -> estimates in mas 

    pdKepEle_FRP(1,1) = dot(pdSatPosGCRF',drdpar_FRP(1,:)');                    % [m/m]     -> [] estimates in cm
    pdKepEle_FRP(2,1) = dot(pdSatPosGCRF',drdpar_FRP(2,:)') *100;               % [cm]      -> estimates dimensionless
    pdKepEle_FRP(3,1) = dot(pdSatPosGCRF',drdpar_FRP(3,:)') *100*(1/rad2mas);   % [cm/mas]  -> estimates in mas
    pdKepEle_FRP(4,1) = dot(pdSatPosGCRF',drdpar_FRP(4,:)') *100*(1/rad2mas);   % [cm/mas]  -> estimates in mas
    pdKepEle_FRP(5,1) = dot(pdSatPosGCRF',drdpar_FRP(5,:)') *100*(1/rad2mas) ;  % [cm/mas]  -> estimates in mas
    pdKepEle_FRP(6,1) = dot(pdSatPosGCRF',drdpar_FRP(6,:)') *100*(1/rad2mas);   % [cm/mas]  -> estimates in mas     
end
