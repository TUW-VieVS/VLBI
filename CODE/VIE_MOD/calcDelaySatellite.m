% ************************************************************************
%   Description:
%   function to calculate the time delay to satellite and to determine the
%   partial derivatives of delay to station coordinates and the satellite
%   position.
%
%   Reference:
%
%   Input:
%       'sources'                    structure array   sources structure array
%       'iSc'                        (1,1)             index of current scan
%       'scan'                       structure array   scan structral array
%       'crsStation1'               (3,1)              CRS coordiantes of station 1
%       'crsStation2'               (3,1)              CRS coordinates of station 2
%       't2c'                        (3,3)             terrrestrial to celestial matrices   
%
%   Output:
%       'tau'                        (1,1)             time delay in TT (Klioner)
%       'scan'                       structure array   scan structral array
%       'crfScPos'                   (3,1)             satellite position at emission time u0
%       'crfScVel'                   (3,1)             satellite velocity at emission time u0
%       'k1a'                        (3,1)             source vector station 1
%       'k2a'                        (3,1)             source vector station 2
%
%   External calls:
%
%   Coded for VieVS:
%   22 November 2021 by H. Wolf - created as external function of vie_mod
%
%   Revision:
%   17 July 2027 by H. Wolf - revised the implmentation of Klioner; added Duev 
%
% ************************************************************************


function [tau, scan, crfScPos, crfScVel, k1a, k2a] = calcDelaySatellite(sources, iSc, scan, crsStation1, crsStation2, t2c)
    global c
    global gme
    global omega
    LG = 6.969290134 *10^(-10);

    ddtThreshold  = 1e-16;     % thresholds for Duev 
    maxIterations = 25;        % max number of iteratoins for Duev
           
    scan_cur = scan(iSc);
    s_cur = sources.s(scan_cur.iso);
  
    [tRefSecInterpol, u1, refidx] = get_reference_epoch_satellite(scan_cur, s_cur);

    pos_crf = [s_cur.x_crf(refidx), s_cur.y_crf(refidx), s_cur.z_crf(refidx)];
    vel_crf = [s_cur.vx_crf(refidx), s_cur.vy_crf(refidx), s_cur.vz_crf(refidx)];
    vel_trf = [s_cur.vx_trf(refidx), s_cur.vy_trf(refidx), s_cur.vz_trf(refidx)];

    % Get spacecraft position at time of observation (CRF):
    crfScPos = lagint9_ultra_fast(tRefSecInterpol, pos_crf, u1);
    
    if s_cur.flag_v_crf
        crfScVel = lagint9_ultra_fast(tRefSecInterpol, vel_crf, u1);
    end

    if s_cur.flag_v_trf
        trfScVel = lagint9_ultra_fast(tRefSecInterpol, vel_trf, u1);
    end
   
    scan(iSc).trfSat = t2c' * crfScPos;
    scan(iSc).crfSat = 1*crfScPos;
    scan(iSc).v_crfSat = crfScVel;
    scan(iSc).v_trfSat = trfScVel;

    % Get spacecraft position at the time of emission u0 (CRF)
    numberOfIterations = 0;
    ddu0 = 999999;
    u0 = u1 - norm(crsStation1 - crfScPos)/c;
    digits(20)
    while(abs(ddu0) > ddtThreshold) 
        numberOfIterations = numberOfIterations + 1;
        u0_old = u0;

        [dugr1, dDelay_dt] = get_dugr(crsStation1, crfScPos, crfScVel);
  
        f = u1 - u0 - norm(crsStation1 - crfScPos)/c - dugr1;
        df = crfScVel'*(crsStation1 - crfScPos)/(c*norm(crsStation1 - crfScPos)) -1 - dDelay_dt;
        u0 = vpa(u0_old - f/df);

        % Get spacecraft position and velocity at time at current estimated emission time u0
        crfScPos = lagint9_ultra_fast(tRefSecInterpol, pos_crf, u0);
        crfScVel = lagint9_ultra_fast(tRefSecInterpol, vel_crf, u0);
        
        [dugr1, ~] = get_dugr(crsStation1, crfScPos, crfScVel);
    
        ddu0 = u1 - u0 - norm(crsStation1 - crfScPos)/c - dugr1;

        if numberOfIterations >= maxIterations
            fprintf(' Warning: Max. number of iterations (%d) for near field delay reached! ddt = %f ps', maxIterations, ddu0*10^(12));
            break;
        end
    end

    % for checks:
    % res = u1 - u0 - norm(crsStation1 - crfScPos)/c - dugr1;
    % res_range = res*c;

    %% Iteration Light time equation 2 
    numberOfIterations_u2 = 0;
    ddu2 = 999999;
    u2 = u0 + norm(crfScPos - crsStation2)/c;

    crsStation2_u2 = crsStation2;

    stationVel = [-omega*crsStation2_u2(2);
                   omega*crsStation2_u2(1);
                   0];

    while(abs(ddu2) > ddtThreshold)
        numberOfIterations_u2 = numberOfIterations_u2 + 1;
        u2_old = u2;

        [dugr2, dDelay_dt] = get_dugr(crsStation2_u2, crfScPos, crfScVel);

        f = u2 - u0 - norm(crsStation2_u2 - crfScPos)/c - dugr2;
        df = 1 - stationVel'*(crsStation2_u2-crfScPos)/(c*norm(crsStation2_u2-crfScPos)) - dDelay_dt;
        u2 = u2_old - f/df;

        delta_t = u2 - u1;

        % Earth rotation during signal propagation
        Rz = [ cos(omega*delta_t)  -sin(omega*delta_t) 0;
              sin(omega*delta_t)  cos(omega*delta_t) 0;
               0                   0                  1];

        crsStation2_u2 = Rz * crsStation2;

        stationVel = [-omega*crsStation2_u2(2);
                       omega*crsStation2_u2(1);
                       0];

        [dugr2, ~] = get_dugr(crsStation2_u2, crfScPos, crfScVel);
        ddu2 = u2 - u0 - norm(crsStation2_u2 - crfScPos)/c - dugr2;

        if numberOfIterations_u2 >= maxIterations
            fprintf(' Warning: Max. number of iterations (%d) for near field delay reached! ddu2 = %f ps \n', maxIterations, ddu2*10^(12));
            break;
        end
    end
 
    % for checks:
    % res = u2 - u0 - norm(crsStation2_u2 - crfScPos)/c - dugr2;
    % res_range = res*c;

    R1 = norm(crfScPos - crsStation1);
    R2 = norm(crfScPos - crsStation2_u2);

    % Vector station-spacecraft (source vectors) at the time of signal emission (spacecraft) and reception at station one (stations)
    L1  = crfScPos' - crsStation1';
    L2  = crfScPos' - crsStation2';
    
    v2 = [-omega*crsStation2(2);
           omega*crsStation2(1);
           0];

    We = gme/norm(crsStation2);

    du0     = (norm(L2) - norm(L1))/c; % Difference in travel time not considering retarded BL effect [sec]
    n       = L2/norm(L2);
    dugr    = 2*gme/c^3*log(...
                ((norm(crsStation2_u2) + norm(crfScPos) + norm(crsStation2_u2 - crfScPos)) * (norm(crsStation1) + norm(crfScPos) - norm(crsStation1 - crfScPos)))...
                /((norm(crsStation2_u2) + norm(crfScPos) - norm(crsStation2_u2 - crfScPos)) * (norm(crsStation1) + norm(crfScPos) + norm(crsStation1 - crfScPos)))); % Gravitaional effect on travel time [sec]

    du_klioner_tcg  = du0*(1-n*v2/c) + dugr;
    du_klioner_tt0  = du0*(1-n*v2/c) - 1/c^2*((v2'*v2)/2+We) * du0  + dugr ; % in TT

    % Duev
    tau_duev_tcg = u2-u1;
    tau_duev_geom_tcg = (R2 - R1)/c  - dugr1 + dugr2;
    tau_duev_geom_tt = tau_duev_geom_tcg * (1 - LG);
    tau_duev_tt = tau_duev_tcg * (1 - LG);

    % Vacuum delay:
    tau = double(du_klioner_tt0);
  
    % diff_duev = (tau_duev_geom_tt - tau_duev_tt)*1e12;
    diff_duev_klioner = (du_klioner_tt0 - tau_duev_tt)*1e12;

    if diff_duev_klioner > 1
        error(' Difference between Klioner and Duev is more than 1 ps! The exact value is: %.8f ps', diff_duev_klioner);
    end

    % Source vectors:
    k1a = double(L1); % Source vector station 1
    k2a = double(L2); % Source vector station 2
end

function [dugr, dDelay_dt] = get_dugr(crsStation, crfScPos, crfScVel)
    global c
    global gme
    rSta = norm(crsStation);
    rSc  = norm(crfScPos);
    rho  = norm(crsStation - crfScPos);
    
    dr_dt   = dot(crfScPos, crfScVel) / rSc;
    drho_dt = dot(crfScPos - crsStation, crfScVel) / rho;
    
    dDelay_dt = 2*gme/c^3 * ( ...
        (dr_dt + drho_dt)/(rSta + rSc + rho) ...
      - (dr_dt - drho_dt)/(rSta + rSc - rho) );

    dugr    = 2*gme/c^3*log(...
            (rSta + rSc + rho )...
            /((rSta + rSc - rho))); 
end
  