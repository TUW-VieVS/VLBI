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


function [ps1, ps2, pdSatPosRSW, pdSatPosGCRF, pdSatPosTRF, tau] = calcDelaySatellite(sources, iSc, scan, idStation1, idStation2, flag_fixSatPostoStat1, ddtThreshold, crsStation1, crsStation2, antenna, secOfDay, mjd, maxIterations, ephem, v2, t2c)
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
        % Get reference epoch for interpolation of SC pos. + vel.:
        refidx  = find((sources.s(scan(iSc).iso).day == scan(iSc).tim(3)) & (sources.s(scan(iSc).iso).sec_of_day > secOfDay));
        if refidx(1) < 8
            error('S.C. ephemeris data do not cover the required time (earlier epochs needed)! Add missing data to S.C. ephem. file!');
        elseif refidx(1)+6 > length(sources.s(scan(iSc).iso).day)
            error('S.C. ephemeris data do not cover the required time (later epochs needed)! Add missing data to S.C. ephem. file!');
        end
        refidx             = refidx(1)-7 : 1 : refidx(1)+6; % suitable for interpolation with "lagint9.m" and "dt" < 1 sec
        tRefSecInterpol  = sources.s(scan(iSc).iso).sec_of_day(refidx);
        tIntegerMjd       = floor(sources.s(scan(iSc).iso).mjd(refidx));
        tRefMjd           = tIntegerMjd(1);
        offsetSec          = (tIntegerMjd - tRefMjd) * 86400; % full days since first interpolation epoch in [sec]
        tRefSecInterpol  = tRefSecInterpol + offsetSec; % add since first interpolation epoch in [sec]
        tRefSec           = tRefSecInterpol(1);
        tRefSecInterpol  = tRefSecInterpol - tRefSec;

        tIntegerMjdObs   = floor(mjd);
        tRefOffsetObs    = (tIntegerMjdObs - tRefMjd) * 86400;
        tRefSecObs       = secOfDay + tRefOffsetObs;
        tRefSecObs       = tRefSecObs -  tRefSec;

        % Get spacecraft position at time of observation (CRF):
        crfScPosX = lagint9(tRefSecInterpol, sources.s(scan(iSc).iso).x_crf(refidx), tRefSecObs);
        crfScPosY = lagint9(tRefSecInterpol, sources.s(scan(iSc).iso).y_crf(refidx), tRefSecObs);
        crfScPosZ = lagint9(tRefSecInterpol, sources.s(scan(iSc).iso).z_crf(refidx), tRefSecObs);
        crfScPos = [crfScPosX; crfScPosY; crfScPosZ];    % (3,1), [m]   

        % Get spacecraft velocity at time of observation (CRF):
        % - Calculated from positions one second before and after the current obs. epoch
        % => Cont. for 2 sec => no further interpolation done
        if ~sources.s(scan(iSc).iso).flag_v_trf          
            crfScPosXm  = lagint9(tRefSecInterpol, sources.s(scan(iSc).iso).x_crf(refidx), tRefSecObs - 1);
            crfScPosYm  = lagint9(tRefSecInterpol, sources.s(scan(iSc).iso).y_crf(refidx), tRefSecObs - 1);
            crfScPosZm  = lagint9(tRefSecInterpol, sources.s(scan(iSc).iso).z_crf(refidx), tRefSecObs - 1);
            crfScPosXp  = lagint9(tRefSecInterpol, sources.s(scan(iSc).iso).x_crf(refidx), tRefSecObs + 1);
            crfScPosYp  = lagint9(tRefSecInterpol, sources.s(scan(iSc).iso).y_crf(refidx), tRefSecObs + 1);
            crfScPosZp  = lagint9(tRefSecInterpol, sources.s(scan(iSc).iso).z_crf(refidx), tRefSecObs + 1);
            crfScVel  = ([crfScPosXp; crfScPosYp; crfScPosZp] - [crfScPosXm; crfScPosYm; crfScPosZm])/2;   %(3,1), [m/sec]
        end

        % Get spacecraft position at the time of emission (CRF):
        % - According to approach "geocneu = 0" in vie_mod_tie.m (lines 735-763) by L. Plank
        crfScPosTmp(1, :) = crfScPos';
        t_ref_sec_obs_tmp = tRefSecObs;

        % Iteration init.:
        numberOfIterations      = 0;
        dt                      = 999999;
        ddt                     = 999999;

        while(abs(ddt) > ddtThreshold) 
            numberOfIterations = numberOfIterations + 1;

            % Get spacecraft velocity at time of observation (CRF):
            % - from ephem. file, if available
            if sources.s(scan(iSc).iso).flag_v_trf
                crfScVelX = lagint9(tRefSecInterpol, sources.s(scan(iSc).iso).vx_crf(refidx), t_ref_sec_obs_tmp);
                crfScVelY = lagint9(tRefSecInterpol, sources.s(scan(iSc).iso).vy_crf(refidx), t_ref_sec_obs_tmp);
                crfScVelZ = lagint9(tRefSecInterpol, sources.s(scan(iSc).iso).vz_crf(refidx), t_ref_sec_obs_tmp);
                crfScVel = [crfScVelX; crfScVelY; crfScVelZ];
            end

            % Correction:
            dt_old  = dt;
            dt      = norm(crfScPos - crsStation1)/c - ((crfScPos' - crsStation1')*crfScVel)/c^2; % [sec] (6.4)
            t_ref_sec_obs_tmp = tRefSecObs - dt; % corrected epoch [sec] since "ref. time" (t_ref_sec, t_ref_mjd)

            %iteration1
            crfScPosX = lagint9(tRefSecInterpol, sources.s(scan(iSc).iso).x_crf(refidx), t_ref_sec_obs_tmp);
            crfScPosY = lagint9(tRefSecInterpol, sources.s(scan(iSc).iso).y_crf(refidx), t_ref_sec_obs_tmp);
            crfScPosZ = lagint9(tRefSecInterpol, sources.s(scan(iSc).iso).z_crf(refidx), t_ref_sec_obs_tmp);
            crfScPos = [crfScPosX; crfScPosY; crfScPosZ];
            crfScPosTmp(numberOfIterations +1, :) = crfScPos';

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

    %K   = (L1+L2)/(norm(L1)+norm(L2)); % mittlere Richtungsvektor

    % Vacuum delay:
    tau = du;

    % Source vectors:
    %k1a = L1; % Source vector station 1 
    %k2a = L2; % Source vector station 2 

    %rqu = K;
    %rq  = K;

    % #### partial derivative of the delay w.r.t. position of the space craft:  ####
    nL1 = norm(L1);
    nL2 = norm(L2);

    % analytical PD in GCRF:
    dudws_part1 = (crfScPos - crsStation2) ./ nL2   -  (crfScPos - crsStation1) ./ nL1; 
    dudws_part2 = v2 - (crfScPos - crsStation1)./nL1 .* 1/nL2 .* (L2*v2)  +  nL1 .* (crfScPos - crsStation2) * 1/nL2^3 * (L2*v2)  - (nL1 * 1/nL2) .* v2; 
    pdSatPosGCRF = dudws_part1./ c - dudws_part2./ c^2; % [sec/m] % PD of delay time du w.r.t. satellite pos. in GCRF [sec/m], (3x1 vector)
    pdSatPosGCRF = pdSatPosGCRF * c; % * 100 / 100; % Unit conversion: [sec/m] => [cm/cm] = []; => estimates will be in [cm]

    % Rotation of PD to RSW system:
    [~, ~, transmat] = rv2rsw(crfScPos, crfScVel);
    pdSatPosRSW = transmat*pdSatPosGCRF;
    
    % Rotation of PD to TRF system:
    pdSatPosTRF = t2c' * pdSatPosGCRF;
        
    % Nummerical derivation of du0 w.r.t. ws1 (du0/dws1)
    s = 0.5;
    crfScPosp = crfScPos + [s; 0; 0]; % plus s m in ws1
    L1p     = crfScPosp' - crsStation1';
    L2p     = crfScPosp' - crsStation2';
    crfScPosm = crfScPos - [s; 0; 0]; % minus s m in ws1
    L1m     = crfScPosm' - crsStation1';
    L2m     = crfScPosm' - crsStation2';
    dudws1_num     = ( ((norm(L2p) - norm(L1p))) - ((norm(L2m) - norm(L1m))) ) / (2*s);

    crfScPosp = crfScPos + [0; s; 0]; % plus s m in ws2
    L1p     = crfScPosp' - crsStation1';
    L2p     = crfScPosp' - crsStation2';
    crfScPosm = crfScPos - [0; s; 0]; % minus s m in ws2
    L1m     = crfScPosm' - crsStation1';
    L2m     = crfScPosm' - crsStation2';
    dudws2_num     = ( ((norm(L2p) - norm(L1p))) - ((norm(L2m) - norm(L1m))) ) / (2*s);

    crfScPosp = crfScPos + [0; 0; s]; % plus s m in ws3
    L1p     = crfScPosp' - crsStation1';
    L2p     = crfScPosp' - crsStation2';
    crfScPosm = crfScPos - [0; 0; s]; % minus s m in ws3
    L1m     = crfScPosm' - crsStation1';
    L2m     = crfScPosm' - crsStation2';
    dudws3_num     = ( ((norm(L2p) - norm(L1p))) - ((norm(L2m) - norm(L1m))) ) / (2*s);
    
    npdSatPosGCRF = [dudws1_num; dudws2_num; dudws3_num];
    
    % check if analytical and numerical solutions are equal
    gcrDfDiffPd = pdSatPosGCRF - npdSatPosGCRF; 
    if gcrDfDiffPd(1,1)> 0.00001 || gcrDfDiffPd(2,1)> 0.00001 || gcrDfDiffPd(3,1)> 0.00001
       disp('Analytical and numerical derivatives of tau wrt satellite position are not equal!')
    end

    % partial derivative of du0 w.r.t. station coordinates (in TRF!)
    ps1 = -t2c'*(L1'/norm(L1)); % partial derivative wrt to station 1
    ps2 = t2c'*(L2'/norm(L2));  % partial derivative wrt to station 2     

end

