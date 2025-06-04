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
%       'idStation1'                 (1,1)             index of Station 1
%       'idStation2'                 (1,1)             index of Station 2
%       'flag_fixSatPostoStat1'      (1,1)             flag, if Satellite Position should be fixed to Station1
%       'ddtThreshold'               (1,1)             threshold for iteration of dt
%       'crsStation1'               (3,1)              CRS coordiantes of station 1
%       'crsStation2'               (3,1)              CRS coordinates of station 2
%       'antenna'                    structure array   antenna structure array
%       'secOfDay'                   (1,1)             seconds of day of this scan
%       'mjd'                        (1,1)             mjd of current Scan
%       'maxIterations'              (1,1)             maximum number ofiterations
%       'ephem'                      structure array   Ephermerides
%       'v2'                         (3,1)             velocity of station 2
%       't2c'                        (3,3)             terrrestrial to celestial matrices   
%
%   Output:
%       'ps1'                        (3,1)             partial derivative of delay wrt to coord of station1 in TRF
%       'ps2'                        (3,1)             partial derivative of delay wrt to coord of station2 in TRF
%       'pdSatPosRSW'                (3,3)             partial derivatives of delay wrt satellite position in RSW-frame
%       'pdSatPosGCRF'               (3,3)             partial derivatives of delay wrt satellite position in GCRF-frame
%       'pdSatPosTRF'                (3,3)             partial derivatives of delay wrt satellite position in TRF-frame
%       'tau'                        (1,1)             time delay
%
%   External calls:
%
%   Coded for VieVS:
%   22 November 2021 by H. Wolf - created as external function of vie_mod
%
%   Revision:
%
% ************************************************************************


function [ps1, ps2, pdSatPosRSW, pdSatPosNTW, pdSatPosGCRF, pdSatPosTRF, tau, scan, k1a, k2a, pGammaSun, fac1] = calcDelaySatellite(sources, iSc, scan, idStation1, idStation2, flag_fixSatPostoStat1, ddtThreshold, crsStation1, crsStation2, antenna, secOfDay, mjd, maxIterations, ephem, v2, t2c)
    % Init.:
    pGammaSun = []; % Not yet calculated for satellite scans
    global c
    global gmm
    global gme

    % Ephemerides:
    earth   = ephem.earth(iSc).xbar;
    sun     = ephem.sun(iSc).xgeo;
    moon    = ephem.moon(iSc).xgeo;

    if (idStation1 == 1) || ~flag_fixSatPostoStat1                 
        scan_cur = scan(iSc);
        s_cur = sources.s(scan_cur.iso);
        % Get reference epoch for interpolation of SC pos. + vel.:
        cur_date = datetime(scan_cur.tim(1), scan_cur.tim(2), scan_cur.tim(3));
        cur_day = cur_date.Day;
        sec_of_day_vec = s_cur.sec_of_day;
        day_vec = s_cur.day;
        orbit_type = s_cur.orbit_file_type;
        
        refidx = find((day_vec == cur_day) & (sec_of_day_vec >= secOfDay));
        if secOfDay > 86100 && ismember(orbit_type, {'tle', 'sp3', 'fso'})
            [sort_sec, sort_idx] = sort(sec_of_day_vec);
            [~, mid_idx] = sort(sort_idx(end-2:end));
            next_date = cur_date + days(1);              
            refidx = find((day_vec == next_date.Day) & (sec_of_day_vec == sort_sec(mid_idx(2))));
        end
        
        if secOfDay == 0 && ismember(orbit_type, {'sp3', 'fso'})
            [~, idx] = sort(sec_of_day_vec);
            mid_idx = idx(end-1);  % mittlerer der letzten drei
            
            if scan_cur.tim(3) == 1
                year_before = scan_cur.tim(1) - 1;
                is_leap = mod(year_before, 400) == 0 || (mod(year_before, 4) == 0 && mod(year_before, 100) ~= 0);
                last_doy = 366 * is_leap + 365 * (~is_leap);
                date_bef = datetime(year_before, 1, last_doy);
            else
                date_bef = cur_date - days(1);
            end
            refidx = find((day_vec == date_bef.Day) & ((sec_of_day_vec == sec_of_day_vec(mid_idx))));
        end
       
        if refidx(1) < 8
            error('S.C. ephemeris data do not cover the required time (earlier epochs needed)! Add missing data to S.C. ephem. file!');
        elseif refidx(1)+6 > length(s_cur.day)
            error('S.C. ephemeris data do not cover the required time (later epochs needed)! Add missing data to S.C. ephem. file!');
        end
        
        nSamples = length(s_cur.sec_of_day);
        refStart = max(refidx(1)-7, 1);
        refEnd = min(refidx(1)+6, nSamples);
        refidx = refStart:refEnd;
        
        tRefSecInterpol = s_cur.sec_of_day(refidx);
        tIntegerMjd = floor(s_cur.mjd(refidx));
        tRefMjd = tIntegerMjd(1);
        offsetSec = (tIntegerMjd - tRefMjd) * 86400;
        tRefSecInterpol = tRefSecInterpol + offsetSec;
        tRefSec = tRefSecInterpol(1);
        tRefSecInterpol = tRefSecInterpol - tRefSec;
        
        tIntegerMjdObs = floor(mjd);
        tRefOffsetObs = (tIntegerMjdObs - tRefMjd) * 86400;
        tRefSecObs = secOfDay + tRefOffsetObs;
        tRefSecObs = tRefSecObs - tRefSec;
      
        x_crf = s_cur.x_crf(refidx);
        y_crf = s_cur.y_crf(refidx);
        z_crf = s_cur.z_crf(refidx);
        vx_crf = s_cur.vx_crf(refidx);
        vy_crf = s_cur.vy_crf(refidx);
        vz_crf = s_cur.vz_crf(refidx);
        vx_trf = s_cur.vx_trf(refidx);
        vy_trf = s_cur.vy_trf(refidx);
        vz_trf = s_cur.vz_trf(refidx);

        % Get spacecraft position at time of observation (CRF):
        crfScPosX = lagint9(tRefSecInterpol, x_crf, tRefSecObs);
        crfScPosY = lagint9(tRefSecInterpol, y_crf, tRefSecObs);
        crfScPosZ = lagint9(tRefSecInterpol, z_crf, tRefSecObs);
        crfScPos = [crfScPosX; crfScPosY; crfScPosZ];    % (3,1), [m]  
        
        if s_cur.flag_v_crf
            crfScVelX = lagint9(tRefSecInterpol, vx_crf, tRefSecObs);
            crfScVelY = lagint9(tRefSecInterpol, vy_crf, tRefSecObs);
            crfScVelZ = lagint9(tRefSecInterpol, vz_crf, tRefSecObs);
            crfScVel = [crfScVelX; crfScVelY; crfScVelZ];
        end

        if s_cur.flag_v_trf
            trfScVelX = lagint9(tRefSecInterpol, vx_trf, tRefSecObs);
            trfScVelY = lagint9(tRefSecInterpol, vy_trf, tRefSecObs);
            trfScVelZ = lagint9(tRefSecInterpol, vz_trf, tRefSecObs);
            trfScVel = [trfScVelX; trfScVelY; trfScVelZ];
        end

        scan(iSc).trfSat = t2c' * crfScPos;
        scan(iSc).crfSat = 1*crfScPos;
        scan(iSc).v_crfSat = crfScVel;
        scan(iSc).v_trfSat = trfScVel;


        % Get spacecraft position at the time of emission (CRF):
        % - According to approach "geocneu = 0" in vie_mod_tie.m (lines 735-763) by L. Plank

        % Iteration init.:
        %crfScPosTmp(1, :) = crfScPos';
        t_ref_sec_obs_tmp = tRefSecObs;
        numberOfIterations      = 0;
        dt                      = 999999;
        ddt                     = 999999;

        while(abs(ddt) > ddtThreshold) 
            numberOfIterations = numberOfIterations + 1;

            % Get spacecraft velocity at time of observation (CRF)
            if s_cur.flag_v_trf && numberOfIterations>1
                crfScVelX = lagint9(tRefSecInterpol, vx_crf, t_ref_sec_obs_tmp);
                crfScVelY = lagint9(tRefSecInterpol, vy_crf, t_ref_sec_obs_tmp);
                crfScVelZ = lagint9(tRefSecInterpol, vz_crf, t_ref_sec_obs_tmp);
                crfScVel = [crfScVelX; crfScVelY; crfScVelZ];
            end

            % Correction
            dt_old  = dt;
            dt      = norm(crfScPos - crsStation1)/c - ((crfScPos' - crsStation1')*crfScVel)/c^2; % [sec] (6.4)
            t_ref_sec_obs_tmp = tRefSecObs - dt; % corrected epoch [sec] since "ref. time" (t_ref_sec, t_ref_mjd)

            %iteration
            crfScPosX = lagint9(tRefSecInterpol, x_crf, t_ref_sec_obs_tmp);
            crfScPosY = lagint9(tRefSecInterpol, y_crf, t_ref_sec_obs_tmp);
            crfScPosZ = lagint9(tRefSecInterpol, z_crf, t_ref_sec_obs_tmp);
            crfScPos = [crfScPosX; crfScPosY; crfScPosZ];
            %crfScPosTmp(numberOfIterations +1, :) = crfScPos';

            ddt = dt_old - dt;
            if numberOfIterations >= maxIterations
                fprintf(' Warning: Max. number of iterations (%d) for near filed delay reached! ddt = %f sec', maxIterations, ddt);
                break;
            end
        end

    end   

    % Vector station-spacecraft (source vectors) at the time of signal emission (spacecraft) and reception at station one (stations)
    L1  = crfScPos' - crsStation1';  % (1x3)
    L2  = crfScPos' - crsStation2';

    % gravitational potential @geocentre / all except earth
    sunb = sun + earth;
    Wsun  = ephem.gms/norm(earth-sunb);
    moonb = moon + earth;
    Wmoon = gmm/norm(earth-moonb);
    Wmerc = ephem.gmmerc/norm(earth-ephem.merc(iSc).xbar);
    Wvenu = ephem.gmvenu/norm(earth-ephem.venu(iSc).xbar);
    Wmars = ephem.gmmars/norm(earth-ephem.mars(iSc).xbar);
    Wjupi = ephem.gmjupi/norm(earth-ephem.jupi(iSc).xbar);
    Wsatu = ephem.gmsatu/norm(earth-ephem.satu(iSc).xbar);
    Wuran = ephem.gmuran/norm(earth-ephem.uran(iSc).xbar);
    Wnept = ephem.gmnept/norm(earth-ephem.nept(iSc).xbar);
    Wplut = ephem.gmplut/norm(earth-ephem.plut(iSc).xbar);
    We = Wplut + Wnept + Wuran + Wsatu + Wjupi + Wmars + Wvenu + Wmerc + Wmoon + Wsun;

    du0     = (norm(L2) - norm(L1))/c; % Difference in travel time not considering retarded BL effect [sec]
    n       = L2/norm(L2);
    dugr    = 2*gme/c^3*log(...
                ((norm(crsStation2) + norm(crfScPos) + norm(crsStation2 - crfScPos)) * (norm(crsStation1) + norm(crfScPos) - norm(crsStation1 - crfScPos)))...
                /((norm(crsStation2) + norm(crfScPos) - norm(crsStation2 - crfScPos)) * (norm(crsStation1) + norm(crfScPos) + norm(crsStation1 - crfScPos)))); % Gravitaional effect on travel time [sec]

    % For the GEOCENTR dugr becomes -inf!
    if strcmp(antenna(idStation1).name,'GEOCENTR') || strcmp(antenna(idStation2).name,'GEOCENTR')
        dugr = 0;
        du  = du0*(1-n*v2/c);% + dugr - 1/c^2*((v2'*v2)/2+We) * du0; % Delay:  Klioner, 1991, formular (6.3) [sec]
    else
        du  = du0*(1-n*v2/c) + dugr - 1/c^2*((v2'*v2)/2+We) * du0; % Delay:  Klioner, 1991, formular (6.3) [sec]
    end

    K   = (L1+L2)/(norm(L1)+norm(L2)); % mittlere Richtungsvektor

    % Vacuum delay:
    tau = du;

    rq  = K;
    rqu = K/norm(K);

    % #### partial derivative of the delay w.r.t. position of the space craft:  ####
    nL1 = norm(L1);
    nL2 = norm(L2);

    % analytical PD in GCRF:
    dudws_part1 = (crfScPos - crsStation2) ./ nL2   -  (crfScPos - crsStation1) ./ nL1; 
    dudws_part2 = v2 - (crfScPos - crsStation1)./nL1 .* 1/nL2 .* (L2*v2)  +  nL1 .* (crfScPos - crsStation2) * 1/nL2^3 * (L2*v2)  - (nL1 * 1/nL2) .* v2; 
    pdSatPosGCRF = dudws_part1./ c - dudws_part2./ c^2; % [sec/m] % PD of delay time du w.r.t. satellite pos. in GCRF [sec/m], (3x1 vector)
    pdSatPosGCRF = pdSatPosGCRF * c; % * 100 / 100; % Unit conversion: [sec/m] => [cm/cm] = []; => estimates will be in [cm]

    % Rotation of PD to RSW system:
    [~, ~, transmatRSW] = rv2rsw(crfScPos, crfScVel);
    pdSatPosRSW = transmatRSW*pdSatPosGCRF;

    % Rotation of PD to NTW system:
    [~, ~, transmatNTW] = rv2ntw(crfScPos, crfScVel);
    pdSatPosNTW = transmatNTW*pdSatPosGCRF;
    
    % Rotation of PD to TRF system:
    pdSatPosTRF = t2c' * pdSatPosGCRF;
        
    %partial derivative of du0 w.r.t. station coordinates (in TRF!)
    ps1 = +t2c'*(L1'/norm(L1)); % partial derivative wrt to station 1 
    ps2 = -t2c'*(L2'/norm(L2));  % partial derivative wrt to station 2
  
    pGammaSun = 0;
    fac1 = 0;
    
    % Source vectors:
    k1a = L1; % Source vector station 1 : crfScPos' - crsStation1'
    k2a = L2; % Source vector station 2 : crfScPos' - crsStation2'
end