% ************************************************************************
%   Description:
%   function to calculate the time delay to satellite and to determine the
%   partial derivatives of delay to station coordinates and the satellite
%   position.
%
%   Reference:
%
%   Input:
%       'sources'                   structure array   sources structure array
%       'iSc'                       (1,1)             index of current scan
%       'scan'                      structure array   scan structral array
%       'crsStation1'               (3,1)             CRS coordiantes of station 1
%       'crsStation2'               (3,1)             CRS coordinates of station 2
%       't2c'                       (3,3)             terrrestrial to celestial matrices
%       'v2'                        (3,3)             velocity of station 2 (CRS)
%       
%
%   Output:
%       'tau'                        (1,1)             time delay in TT (Klioner)
%       'scan'                       structure array   scan structural array
%       'crfSatPos'                   (3,1)             satellite position at emission time u0
%       'crfSatVel'                   (3,1)             satellite velocity at emission time u0
%       'k1a'                        (3,1)             source vector station 1
%       'k2a'                        (3,1)             source vector station 2
%
%   External calls:
%
%   Coded for VieVS:
%   22 November 2021 by H. Wolf - created as external function of vie_mod
%
%   Revision:
%   17 July 2026 by H. Wolf - revised the implementation of Klioner; added Duev 
%   11 Nov  2026 by H. Wolf - revised the code; fixed TCG and TT inconsistency
%
% ************************************************************************


function [tau, scan, crfSatPos, crfSatVel, k1a, k2a] = calc_nfd(sources, iSc, scan, crsStation1, crsStation2, t2c, v2)
    global c omega
    gme_tt = 3.986004415e14; % [m^3/s^2], TT-compatible
    omegaCrf = t2c * [0; 0; omega];

    RESIDUAL_TOL = 1e-14;      % [s]  convergence threshold of light-time equations
    MAX_ITER     = 25;         % max number of iterations for lighttime equation of Duev
    MAX_DIFF_PS  = 1;          % [ps] allowed deviation analytic vs. iterative
           
    curScan = scan(iSc);
    sat = sources.s(curScan.iso);
    [tRefInterp, u1, refIdx] = get_reference_epoch_satellite(curScan, sat);

    posCrf = [sat.x_crf(refIdx), sat.y_crf(refIdx), sat.z_crf(refIdx)];
    velCrf = [sat.vx_crf(refIdx), sat.vy_crf(refIdx), sat.vz_crf(refIdx)];
    velTrf = [sat.vx_trf(refIdx), sat.vy_trf(refIdx), sat.vz_trf(refIdx)];

    % Get spacecraft position at time of observation (CRF):
    crfSatPos = lagint9_fast(tRefInterp, posCrf, u1);
    
    if sat.flag_v_crf
        crfSatVel = lagint9_fast(tRefInterp, velCrf, u1);
    end
    if sat.flag_v_trf
        trfSatVel = lagint9_fast(tRefInterp, velTrf, u1);
    end
   
    scan(iSc).trfSat = t2c' * crfSatPos;
    scan(iSc).crfSat = crfSatPos;
    scan(iSc).v_crfSat = crfSatVel;
    scan(iSc).v_trfSat = trfSatVel;

    %% Light-time equation 1: lightTime1 = u1 - u0 
    % Get spacecraft position at the time of emission u0 (station1 fixed)
    lightTime1   = norm(crsStation1 - crfSatPos)/c;
    
    for iter1 = 1:MAX_ITER
        u0     = u1 - lightTime1;
        crfSatPos = lagint9_fast(tRefInterp, posCrf, u0);
        crfSatVel = lagint9_fast(tRefInterp, velCrf, u0);

        rangeVec = crsStation1 - crfSatPos;
        rho      = norm(rangeVec);
        [shapiro1, dShapiroDt] = get_shapiro_delay(crsStation1, crfSatPos, crfSatVel, gme_tt);

        residual1 = lightTime1 - rho/c - shapiro1;
        if abs(residual1) <= RESIDUAL_TOL, break; end

        fDeriv     = crfSatVel' * rangeVec / (c*rho) - 1 - dShapiroDt;
        lightTime1 = lightTime1 + residual1 / fDeriv;

    end
    if abs(residual1) > RESIDUAL_TOL
        warning('calc_nfd:noConvergence', 'Light Time Equation 1: %d iterations reached, residual = %.3f ps', MAX_ITER, residual1*1e12);
    end

    %% Iteration Light time equation 2 
    lightTime2  = norm(crfSatPos - crsStation2) / c;
    for iter2 = 1:MAX_ITER
        rotAngle = omega * (lightTime2 - lightTime1);  % Earth rotation u1 -> u2
        Rz = [cos(rotAngle) -sin(rotAngle) 0;
              sin(rotAngle)  cos(rotAngle) 0;
              0              0             1];
        crsStation2_u2 = t2c * (Rz * (t2c' * crsStation2));
        stationVel     = cross(omegaCrf, crsStation2_u2);

        rangeVec = crsStation2_u2 - crfSatPos;
        rho      = norm(rangeVec);
        [shapiro2, ~] = get_shapiro_delay(crsStation2_u2, crfSatPos, crfSatVel, gme_tt);

        residual2 = lightTime2 - rho/c - shapiro2;
        if abs(residual2) <= RESIDUAL_TOL, break; end

        fDeriv     = 1 - stationVel' * rangeVec / (c*rho); 
        lightTime2 = lightTime2 - residual2 / fDeriv;
    end
    if abs(residual2) > RESIDUAL_TOL
        warning('calc_nfd:noConvergence', 'Light Time Equation 2: %d iterations reached, residual = %.3f ps', MAX_ITER, residual2*1e12);
    end
     
    %% Analytic formula (Klioner)
    L1  = crfSatPos' - crsStation1';
    L2  = crfSatPos' - crsStation2';
    du0     = (norm(L2) - norm(L1))/c; 
    dugr    = 2*gme_tt/c^3*log(...
                ((norm(crsStation2_u2) + norm(crfSatPos) + norm(crsStation2_u2 - crfSatPos)) * (norm(crsStation1) + norm(crfSatPos) - norm(crsStation1 - crfSatPos)))...
                /((norm(crsStation2_u2) + norm(crfSatPos) - norm(crsStation2_u2 - crfSatPos)) * (norm(crsStation1) + norm(crfSatPos) + norm(crsStation1 - crfSatPos)))); % Gravitational effect on travel time [sec]
    tauKlioner = du0 * (1 - dot(L2/norm(L2), v2)/c) + dugr; % [sec] in TT

    %% Duev
    tauDuev = lightTime2-lightTime1; % [sec] in TT
    
    %% Consistency check
    diffPs = (tauKlioner - tauDuev)*1e12;
    if abs(diffPs) > MAX_DIFF_PS
         warning('calcNfd:modelMismatch', ...
        'Scan %d: analytic (Klioner) and iterative (Duev) delay differ by %.8f ps (limit %g ps).', ...
        iSc, diffPs, MAX_DIFF_PS);
    end

    tau = tauKlioner;
    k1a = crfSatPos' - crsStation1';    % Source vector station 1
    k2a = crfSatPos' - crsStation2_u2'; % Source vector station 2
end

function [shapiro, dShapiroDt] = get_shapiro_delay(crsStation, crfSatPos, crfSatVel, gme_tt)
    global c
    rSta = norm(crsStation);
    rSat  = norm(crfSatPos);
    rho  = norm(crsStation - crfSatPos);
    
    drSat   = dot(crfSatPos, crfSatVel) / rSat;
    dRho = dot(crfSatPos - crsStation, crfSatVel) / rho;
    
    dShapiroDt = 2*gme_tt/c^3 * ( ...
        (drSat + dRho)/(rSta + rSat + rho) ...
      - (drSat - dRho)/(rSta + rSat - rho) );

    shapiro    = 2*gme_tt/c^3*log(...
            (rSta + rSat + rho )...
            /((rSta + rSat - rho))); 
end
  