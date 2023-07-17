% ************************************************************************
%   Description:
%   function to correct antenna position due to solid earth tides, tidal
%   ocean loading, tidal atmosphere loading, non-tidal atmosphere loading,
%   GIA uplift rates, hydrology loading, pole tide, ocean pole tide, 
%   seasonal variations. Further the tropospheric paramters are determined
%   and stored in the scan structure.
%
%   Reference:
%
%   Input:
%       'iSc'          (1,1)             index of current scan 
%       'iStat'        (1,1)             index of current station
%       'scan'         structure array   scan structure 
%       'antenna'      structure array   antenna structure 
%       'parameter'    structure array   parameter structure 
%       'session'      string            session name
%       'mjd'          (1,1)             mjd of current Scan            
%       'ANT'          structure array   antenna positions
%       'VEL'          structure array   antenna velocities 
%       't2c'          (3,3)             terrrestrial to celestial matrices
%       'leap'         (3,1)             velocity of station 2
%       'cto_F'        (1,x)             frequencies of tidal constituents 
%       'cto_P'        (1,x)             phases of tidal constituents
%       'cto_TAMP'     (1,x)             Cartwright-Edden amplitudes of tidal constituents
%       'cto_IDD1'     (x,1)             first digit of the Doodson number (used in libiers_hardisp/libiers_admint) 
%       'PHI'          (n,1)             elliptical coordinate of antenna (latitude)
%       'LAM'          (n,1)             elliptical coordinate of antenna (longitude)
%       'tim'          (7,1)             time of current scan 
%       'xp'           (1,1)             x-pole coordinates at time of current scan
%       'yp'           (1,1)             y-pole coordinates at time of current scan
%       'cpsd_all'     (3,y,n)           matrix containing post-seismic deformations for all stations in antenna file for all mjds
%       'ephem'        structure array   Ephermerides
%
%   Output:
%       'scan'           structure array   scan structure 
%       'flgm_ctp'       (1,1)             flagmessage if a linear model instead of a cubic 
%                                          one for the mean pole had to be used             
%
%   External calls:
%        mathews.m      matthew.m        libiers_hardisp.m  ren2xyz.m  
%        ctatmos12.m    call_spline_4.m  ctpole.m           ctoceanpole.m
%        stseasonal.m   load_trpfile.m   get_trpdel.m 
%
%   Coded for VieVS:
%   22 November 2021 by H. Wolf - created as external function of vie_mod
%
%   Revision:
%
% ************************************************************************


function [scan, flgm_ctp] = antennaCorrections(iSc, iStat, scan, antenna, opt, parameter, session, mjd, ANT, VEL, t2c, leap, cto_F, cto_P, cto_TAMP, cto_IDD1, PHI, LAM, tim, xp, yp, cpsd_all, ephem, sourceNames)
    moon    = ephem.moon(iSc).xgeo;
    sun     = ephem.sun(iSc).xgeo;
    % Choose only active stations in the scan:

    flgm_ctp    = [];
    if ~isempty(scan(iSc).stat(iStat).temp)
        raytr_used = strcmp(parameter.vie_init.tropSource.name,'raytr');
        ant = ANT(iStat,:);
        vel = VEL(iStat,:);

        % Solid Earth Tide corrections
        cts = [0 0 0];

        pdxyz_h = [0 0 0];
        pdxyz_l = [0 0 0];
        pdxyz_FCN=[0 0 0];
        if parameter.vie_mod.cts == 1
            if opt.est_love || opt.est_shida || opt.est_FCNset
                [cts, pdxyz_h, pdxyz_l, pdxyz_FCN] = matthew(mjd, leap, t2c, ant, moon, sun);
            else
                [cts] = mathews(mjd, leap, t2c, ant, moon, sun);
            end
        end

        % Tidal Ocean Loading corrections
        cto = [0 0 0];
        if parameter.vie_mod.cto == 1
            if isempty(antenna(iStat).cto) == 0
                cto_data = antenna(iStat).cto;
                %[cto] = ctocean(mjd,leap,ant(ist,:),cto_data);
                [cto_rsw] = libiers_hardisp(mjd, leap, cto_data, cto_F,cto_P, cto_TAMP, cto_IDD1); % [Radial, South, West];  IERS Conv. 2010
                cto = ren2xyz([cto_rsw(1) -cto_rsw(3) -cto_rsw(2)] , PHI(iStat), LAM(iStat));
            end
        end

        % Tidal Atmosphere Loading
        cta = [0 0 0];
        if parameter.vie_mod.cta == 1
            if isempty(antenna(iStat).cta) == 0
                S12_data = antenna(iStat).cta;
                S12_data = S12_data(1:12);         % for "ATIDE_LEONID.mat"
                [cta] = ctatmos12(mjd, ant, S12_data);
            end
        end

        % Non-tidal Atmosphere Loading
        cnta = [0 0 0];
        if parameter.vie_mod.cnta == 1
            if isempty(antenna(iStat).cnta_dx) == 0
                cnta_data = antenna(iStat).cnta_dx;
                cnta(:,1) = call_spline_4(cnta_data(:,1), cnta_data(:,2), mjd);
                cnta(:,2) = call_spline_4(cnta_data(:,1), cnta_data(:,3), mjd);
                cnta(:,3) = call_spline_4(cnta_data(:,1), cnta_data(:,4), mjd);
            end
        end

        crg  = [0 0 0]; %xyz [m]
        cgia = [0 0 0]; %xyz [m]

        % Atmosphere loading from regression coefficients
        pres = scan(iSc).stat(iStat).pres;  % total pressure [hPa]
        pres0 = antenna(iStat).crg(1); %[hPa]
        if parameter.vie_mod.crg == 1
            rg = antenna(iStat).crg(2); %[m/hPa]
            crg_h = (pres - pres0) * rg;
            crg = ren2xyz([crg_h 0 0], PHI(iStat), LAM(iStat));
        end
        pcrg_h      = (pres-pres0); %[hPa]
        pcrg_xyz    = ren2xyz([pcrg_h 0 0],PHI(iStat),LAM(iStat));  %[hPa]

        % GIA uplift rates
        if parameter.vie_mod.gia == 1
            vgia = antenna(iStat).gia_dvx(2:4); % m/y
            mjd0_gia = antenna(iStat).gia_dvx(1); %mjd

            cgia = vgia /365.25 * (mjd - mjd0_gia); %m
        end        

        % Hydrology loading
        chl = [0 0 0];
        if parameter.vie_mod.chl == 1
            if isempty(antenna(iStat).chl_dx) == 0
                chl_data = antenna(iStat).chl_dx;
                chl(:,1) = call_spline_4(chl_data(:,1), chl_data(:,2), mjd);
                chl(:,2) = call_spline_4(chl_data(:,1), chl_data(:,3), mjd);
                chl(:,3) = call_spline_4(chl_data(:,1), chl_data(:,4), mjd);
            end
        end
        
        % Non-tidal ocean loading
        cntol = [0 0 0];
        if parameter.vie_mod.cntol == 1
            if isempty(antenna(iStat).cntol_dx) == 0
                cntol_data = antenna(iStat).cntol_dx;
                cntol(:,1) = call_spline_4(cntol_data(:,1), cntol_data(:,2), mjd);
                cntol(:,2) = call_spline_4(cntol_data(:,1), cntol_data(:,3), mjd);
                cntol(:,3) = call_spline_4(cntol_data(:,1), cntol_data(:,4), mjd);
            end
        end


        % Pole Tide (mean: linear or cubic)
        ctp         = [0 0 0];
        flgm_ctp    = [];
        phpole      = [0 0 0];
        plpole      = [0 0 0];
        if parameter.vie_mod.ctp == 1
            ctpm = parameter.vie_mod.ctpm; % mean: linear or cubic
            [ctp, flgm_ctp, phpole, plpole] = ctpole(tim, ant, xp, yp, ctpm);
        end

        % Ocean Pole Tide (mean: linear or cubic)
        ctop        = [0 0 0];
        flgm_ctp    = [];
        if parameter.vie_mod.ctop == 1
            if isempty(antenna(iStat).opl) == 0
                opl = antenna(iStat).opl;
                ctpm = parameter.vie_mod.ctpm; % mean: linear or cubic
                [ctop, flgm_ctp] = ctoceanpole(tim, ant, xp, yp, opl, ctpm);
            end
        end

        % Seasonal variations of station positions
        % (zero a priori; only partial derivatives w.r.t. amplitudes are build up)
        pAcr_xyz = [0 0 0]; pAce_xyz = [0 0 0]; pAcn_xyz = [0 0 0];
        pAsr_xyz = [0 0 0]; pAse_xyz = [0 0 0]; pAsn_xyz = [0 0 0];
        if opt.est_stsespos
            [pAcr_xyz, pAce_xyz, pAcn_xyz, pAsr_xyz, pAse_xyz, pAsn_xyz]= stseasonal(mjd,PHI(iStat),LAM(iStat));
        end

        % post-seismic deformation (ITRF2014 at least)
        cpsd = cpsd_all(:, iSc,iStat)'; % [3 x nScans x nStat] matrix / meters!

        % Displacement of reference markers on the crust
        % => For the GEOCENTR station displacements are set to zero!
        if strcmp(antenna(iStat).name,'GEOCENTR')
            cposit = 0;
        else
            cposit = cts + cto + cta + cnta + crg + cgia + ctp + ctop + chl + cntol + cpsd;
        end

        % time since reference epoch [years]
        tep = (mjd - antenna(iStat).epoch)/365.25;

        % Real antenna position
        antv     = ant + tep * vel;
        ant_trs  = antv + cposit + antenna(iStat).c_ecc;
        ant_crs  = t2c * ant_trs'; % transformation to the cel. system

        % store results per scan & station
        scan(iSc).stat(iStat).x      = ant_trs; % corrected station position in TRS [x,y,z]
        scan(iSc).stat(iStat).xcrs   = ant_crs; % corrected station position in CRS [x,y,z]

        scan(iSc).stat(iStat).pAcr_xyz = pAcr_xyz; % [-] partials for cosine Ampl. - seasonal stat. position variations
        scan(iSc).stat(iStat).pAce_xyz = pAce_xyz; % [-]
        scan(iSc).stat(iStat).pAcn_xyz = pAcn_xyz; % [-]
        scan(iSc).stat(iStat).pAsr_xyz = pAsr_xyz; % [-] partials for sine Ampl. - seasonal stat. position variations
        scan(iSc).stat(iStat).pAse_xyz = pAse_xyz; % [-]
        scan(iSc).stat(iStat).pAsn_xyz = pAsn_xyz; % [-]

        scan(iSc).stat(iStat).phpole = phpole;     % [cm] partials for pole tide Love numbers
        scan(iSc).stat(iStat).plpole = plpole;     % [cm]

        scan(iSc).stat(iStat).prg = pcrg_xyz;      % [hPa] % partials for APL RG

        scan(iSc).stat(iStat).pLove    = pdxyz_h;  % [cm] partials for Love numbers
        scan(iSc).stat(iStat).pShida   = pdxyz_l;  % [cm] partials for Shida numbers
        scan(iSc).stat(iStat).pFCN     = pdxyz_FCN;% [cm] partials for FCN

        % --------------------------------------------------
        %  tropospheric parameters (using ray-traced delays)
        % --------------------------------------------------

        if raytr_used   % if ray-tracing files are used
            if ~exist('firstRaytrRun', 'var')   % read the .trp-file only in the first run
                [raytrdata, raytrFileFoundLog] = load_trpfile(parameter, session);
            end
            scan = get_trpdel(raytrdata, scan, iSc, iStat, antenna, raytrFileFoundLog, sourceNames);
        end

        aht = 0; awt = 0; zhdt = 0; zwdt = 0;
        if strcmpi(parameter.vie_init.tropSource.name,'indModeling')   ||   isempty(scan(iSc).stat(iStat).trop)

            if strcmpi(parameter.vie_init.zhd,'vmf3')
                tmjd = antenna(iStat).vmf3(:,1);
                zhd  = antenna(iStat).vmf3(:,4);

                zhdt = call_spline_4(tmjd,zhd,mjd);
            end
            if strcmpi(parameter.vie_init.zhd,'vmf1')
                tmjd = antenna(iStat).vmf1(:,1);
                zhd  = antenna(iStat).vmf1(:,4);

                zhdt = call_spline_4(tmjd,zhd,mjd);
            end

            if strcmpi(parameter.vie_init.zwd,'vmf3')
                tmjd = antenna(iStat).vmf3(:,1);
                zwd  = antenna(iStat).vmf3(:,5);

                zwdt = call_spline_4(tmjd,zwd,mjd);
            end
            if strcmpi(parameter.vie_init.zwd,'vmf1')
                tmjd = antenna(iStat).vmf1(:,1);
                zwd  = antenna(iStat).vmf1(:,5);

                zwdt = call_spline_4(tmjd,zwd,mjd);
            end

            if strcmpi(parameter.vie_mod.mfh,'vmf3')
                if isempty(antenna(iStat).vmf3) == 0
                    tmjd = antenna(iStat).vmf3(:,1);
                    ah   = antenna(iStat).vmf3(:,2);
                    aw   = antenna(iStat).vmf3(:,3);

                    aht = call_spline_4(tmjd,ah,mjd);
                    awt = call_spline_4(tmjd,aw,mjd);
                end
            end
            if strcmpi(parameter.vie_mod.mfh,'vmf1')
                if isempty(antenna(iStat).vmf1) == 0
                    tmjd = antenna(iStat).vmf1(:,1);
                    ah   = antenna(iStat).vmf1(:,2);
                    aw   = antenna(iStat).vmf1(:,3);

                    aht = call_spline_4(tmjd,ah,mjd);
                    awt = call_spline_4(tmjd,aw,mjd);
                end
            end

            if strcmpi(parameter.vie_mod.mfw,'vmf3')
                if isempty(antenna(iStat).vmf3) == 0
                    tmjd = antenna(iStat).vmf3(:,1);
                    ah   = antenna(iStat).vmf3(:,2);
                    aw   = antenna(iStat).vmf3(:,3);

                    aht = call_spline_4(tmjd,ah,mjd);
                    awt = call_spline_4(tmjd,aw,mjd);
                end
            end
            if strcmpi(parameter.vie_mod.mfw,'vmf1')
                if isempty(antenna(iStat).vmf1) == 0
                    tmjd = antenna(iStat).vmf1(:,1);
                    ah   = antenna(iStat).vmf1(:,2);
                    aw   = antenna(iStat).vmf1(:,3);

                    aht = call_spline_4(tmjd,ah,mjd);
                    awt = call_spline_4(tmjd,aw,mjd);
                end
            end     
        end % use VMF/GPT3/GMF intead of trp files

        % store results per scan & station
        scan(iSc).stat(iStat).aht  = aht;    % interpolated ah
        scan(iSc).stat(iStat).awt  = awt;    % interpolated aw
        scan(iSc).stat(iStat).zhdt = zhdt;   % interpolated zhd
        scan(iSc).stat(iStat).zwdt = zwdt;   % interpolated zwd
    end
end

