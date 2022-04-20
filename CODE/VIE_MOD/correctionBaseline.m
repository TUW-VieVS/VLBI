% ************************************************************************
%   Description:
%   function to correct the time delay tau due to tropospheric delay, axis
%   offset correction, thermal deformation, gravitational deformation and 
%   to calculate the northern and eastern gradients. Further, the external
%   ionospheric delay is determined and added, if it is set to be done. 
%
%   Reference:
%
%   Input:
%       'scan'             structure array   scan structure            
%       'antenna'          structure array   antenna structure 
%       'parameter'        structure array   parameter structure
%       't2c'              (3,3)             terrrestrial to celestial matrice of current scan
%       'mjd'              (1,1)             mjd of current Scan
%       'iSc'              (1,1)             index of current scan
%       'idStation1'       (1,1)             index of Station 1
%       'idStation2'       (1,1)             index of Station 2
%       'k1a'              (1,3)             aberrated source vector
%       'k2a'              (1,3)             aberrated source vector  
%       'rqu'              (1,3)             unit source vector barycentrum-source
%       'v2'               (3,1)             velocity of station 2
%       'v1'               (3,1)             velocity of station 1
%       'tau'              (1,1)             time delay 
%       'cell_grid_GPT3'   cell grid         gpt3 grid    

%
%   Output:
%       'a_ngr'            (n,1)             north  gradients
%       'a_egr'            (n,1)             east gradients
%       'scan'             structure array   scan structure array
%       'antenna'          structure array   antenna structure array
%       'tau'              (1,1)             time delay
%
%   External calls:
%        xyz2ell.m       locsource.m    axis_stat.m    asknewet.m       vmf3.m 
%        vmf1.m          vmf3_ht.m      gpt3_5_fast.m  call_spline_4.m  
%        load_ionfile.m  get_iondel.m   ao_altcorr.m   dao.m           
%
%   Coded for VieVS:
%   22 November 2021 by H. Wolf - created as external function of vie_mod
%
%   Revision:
%
% ************************************************************************

function [a_ngr, a_egr, scan, antenna, tau] = correctionBaseline(scan, antenna, parameter, t2c, mjd, iSc, idStation1, idStation2, k1a, k2a, rqu, v2, v1, tau, cell_grid_GPT3, session, iobs)

    global c

    % Preallocate a priori gradients for each antenna (north and east gradient)
    a_ngr_h = zeros(length(antenna), 1);
    a_egr_h = zeros(length(antenna), 1);
    a_ngr_w = zeros(length(antenna), 1);
    a_egr_w = zeros(length(antenna), 1);
    a_ngr = zeros(length(antenna), 1);
    a_egr = zeros(length(antenna), 1);

    for kstat=1:2
        switch kstat
            case 1
                stnum = idStation1; ka = k1a; % stt = stt1;
            case 2
                stnum = idStation2; ka = k2a; % stt = stt2;
        end

        % only do once for one station for one scan
        if isempty(scan(iSc).stat(stnum).zd)

            % Get corrected station coordinates for actual observation epoch:
            [phi, lam, hell] = xyz2ell(scan(iSc).stat(stnum).x);
            [azim,zd,~,corz,de,LHAe] = locsource(lam,phi,ka,t2c);

            aname  = antenna(stnum).name;

            % AXIS OFFSET
            axtyp  = antenna(stnum).axtyp;
            offs   = antenna(stnum).offs;
            if ~isempty(axtyp)
                [axkt,daxkt]  = axis_stat(phi,azim,zd,corz,de,axtyp,offs,aname);
            else
                axkt=0;
                daxkt=0;
                flagmess.axtyp(stnum) = 1;
            end

            if strcmpi(parameter.vie_init.tropSource.name,'indModeling')  ||  isempty(scan(iSc).stat(stnum).trop)

                % TROPOSPHERE    
                % apriori zenith delay Lz [m] Marini tropospheric model

                % hydrostatic
                if strcmpi(parameter.vie_init.zhd,'in situ')
                    pres = scan(iSc).stat(stnum).pres;  % [hPa]
                    zdry = 0.0022768*pres / (1-0.00266*cos(2*phi)-(0.28e-6*hell));   %[m]
                elseif strcmpi(parameter.vie_init.zhd,'no')
                    zdry = 0;
                elseif strcmpi(parameter.vie_init.zhd,'vmf3')
                    zdry = scan(iSc).stat(stnum).zhdt;
                elseif strcmpi(parameter.vie_init.zhd,'vmf1')
                    zdry = scan(iSc).stat(stnum).zhdt;
                elseif strcmpi(parameter.vie_init.zhd,'gpt3')
                    pres = antenna(stnum).gpt3.p;
                    zdry = 0.0022768*pres / (1-0.00266*cos(2*phi)-(0.28e-6*hell));   %[m]
                else
                    error('Something is wrong here...');
                end

                % wet
                if strcmpi(parameter.vie_init.zwd,'no')
                    zwet = 0;
                elseif strcmpi(parameter.vie_init.zwd,'in situ')
                    e      = scan(iSc).stat(stnum).e;  % [hPa]
                    Tm     = antenna(stnum).gpt3.Tm;
                    lambda = antenna(stnum).gpt3.lambda;
                    zwet   = asknewet ( e , Tm , lambda );
                elseif strcmpi(parameter.vie_init.zwd,'vmf3')
                    zwet   = scan(iSc).stat(stnum).zwdt;
                elseif strcmpi(parameter.vie_init.zwd,'vmf1')
                    zwet   = scan(iSc).stat(stnum).zwdt;
                elseif strcmpi(parameter.vie_init.zwd,'gpt3')
                    e      = antenna(stnum).gpt3.e;
                    Tm     = antenna(stnum).gpt3.Tm;
                    lambda = antenna(stnum).gpt3.lambda;
                    zwet   = asknewet( e , Tm , lambda );
                else
                    error('Something is wrong here...');
                end       

                % mapping function

                % hydrostatic
                aht = scan(iSc).stat(stnum).aht;
                awt = scan(iSc).stat(stnum).awt;
                if strcmpi(parameter.vie_mod.mfh,'vmf3') == 1 && aht ~= 0
                    [mfh,~] = vmf3(aht,awt,mjd,phi,lam,zd);   % VMF3
                elseif strcmpi(parameter.vie_mod.mfh,'vmf1') == 1 && aht ~= 0
                    [mfh,~] = vmf1(aht,awt,mjd,phi,zd);       % VMF1
                elseif strcmpi(parameter.vie_mod.mfh,'gpt3')
                    [~,~,~,~,~,aht,awt,~,~,~,~,~,~] = gpt3_5_fast (mjd,phi,lam,hell,0,cell_grid_GPT3);
                    [mfh,~] = vmf3_ht(aht,awt,mjd,phi,lam,hell,zd);   % GPT3
                else
                    [mfh,~] = vmf3_ht(antenna(stnum).gpt3.ah,antenna(stnum).gpt3.aw,mjd,phi,lam,hell,zd);   % GPT3 as backup if aht = 0
                end

                % wet
                aht = scan(iSc).stat(stnum).aht;
                awt = scan(iSc).stat(stnum).awt;
                if strcmpi(parameter.vie_mod.mfw,'vmf3') == 1 && awt ~= 0
                    [~,mfw] = vmf3(aht,awt,mjd,phi,lam,zd);   % VMF3
                elseif strcmpi(parameter.vie_mod.mfw,'vmf1') == 1 && awt ~= 0
                    [~,mfw] = vmf1(aht,awt,mjd,phi,zd);   % VMF1
                elseif strcmpi(parameter.vie_mod.mfw,'gpt3')
                    [~,~,~,~,~,aht,awt,~,~,~,~,~,~] = gpt3_5_fast (mjd,phi,lam,hell,0,cell_grid_GPT3);
                    [~,mfw] = vmf3_ht (aht,awt,mjd,phi,lam,hell,zd);   % GPT3
                else
                    [~,mfw] = vmf3_ht (antenna(stnum).gpt3.ah,antenna(stnum).gpt3.aw,mjd,phi,lam,hell,zd);   % GPT3 as backup if aht = 0
                end

                % a priori gradient delay model

                % hydrostatic
                switch parameter.vie_mod.apgm_h
                    case 'no'
                        aprgrd_h = 0;
                    case 'grad'
                        %get only wanted lines (i.e. lines with only session-included stations)
                        wantedLines = ismember(gradData{1} , aname);
                        if min(abs(mjd-gradData{2}(wantedLines))) < 0.5   % only use GRAD if there is data at least one half day before or after
                            a_ngr_h(stnum) = call_spline_4(gradData{2}(wantedLines), gradData{3}(wantedLines), mjd)/1000;   % [m]
                            a_egr_h(stnum) = call_spline_4(gradData{2}(wantedLines), gradData{4}(wantedLines), mjd)/1000;   % [m]
                            aprgrd_h = 1/(sin(pi/2-zd)*tan(pi/2-zd)+0.0031)*(a_ngr_h(stnum)*cos(azim)+a_egr_h(stnum)*sin(azim));   % [m]
                        else
                            aprgrd_h = 0;
                            if antenna(stnum).noGrad ~= 1   % write to command window, but only for first observation of session (flag is valid for one whole session)
                                antenna(stnum).noGrad = 1;
                                fprintf('%s%s%s\n','GRAD data not available for station ',strtrim(aname),' at some epochs of this session, no a priori gradients used instead');
                            end
                        end        
                    case 'gpt3'
                        [~,~,~,~,~,~,~,~,~,Gn_h,Ge_h,~,~] = gpt3_5_fast (mjd,phi,lam,hell,0,cell_grid_GPT3);
                        a_ngr_h(stnum) = Gn_h;   % [m]
                        a_egr_h(stnum) = Ge_h;   % [m]
                        aprgrd_h = 1/(sin(pi/2-zd)*tan(pi/2-zd)+0.0031)*(a_ngr_h(stnum)*cos(azim)+a_egr_h(stnum)*sin(azim));   % [m]
                    case 'dao'
                        [aprgrd_h , a_ngr_h(stnum) , a_egr_h(stnum) ] = dao(aname,azim,pi/2-zd);   % [[m],[mm],[mm]], total gradients, despite being named apgrd_h (MacMillan, 1995)
                        aprgrd_w = 0;
                        a_ngr_h(stnum) = a_ngr_h(stnum)/1000;   % [m]
                        a_egr_h(stnum) = a_egr_h(stnum)/1000;   % [m]
                        a_ngr_w(stnum) = 0;
                        a_egr_w(stnum) = 0;
                        if aprgrd_h==0
                            flagmess.dao(stnum) = 1;
                        end 
                end

                % wet
                if ~strcmpi(parameter.vie_mod.apgm_h,'dao')   % only determine a wet part if other model than DAO is used
                    switch parameter.vie_mod.apgm_w
                        case 'no'
                            aprgrd_w = 0;
                        case 'grad'
                           %get only wanted lines (i.e. lines with only session-included stations)
                            wantedLines = ismember(gradData{1} , aname);
                            if min(abs(mjd-gradData{2}(wantedLines))) < 0.5   % only use GRAD if there is data at least one half day before or after from which can be extrapolated
                                a_ngr_w(stnum) = call_spline_4(gradData{2}(wantedLines), gradData{5}(wantedLines), mjd)/1000;   % [m]
                                a_egr_w(stnum) = call_spline_4(gradData{2}(wantedLines), gradData{6}(wantedLines), mjd)/1000;   % [m]
                                aprgrd_w = 1/(sin(pi/2-zd)*tan(pi/2-zd)+0.0007)*(a_ngr_w(stnum)*cos(azim)+a_egr_w(stnum)*sin(azim));   % [m]
                            else
                                aprgrd_w = 0;
                                if antenna(stnum).noGrad ~= 1   % write to command window, but only for first observation of session (flag is valid for one whole session)
                                    antenna(stnum).noGrad = 1;
                                    fprintf('%s%s%s\n','GRAD data not available for station ',strtrim(aname),' at some epochs of this session, no a priori gradients used instead');
                                end
                            end
                        case 'gpt3'
                            [~,~,~,~,~,~,~,~,~,~,~,Gn_w,Ge_w] = gpt3_5_fast (mjd,phi,lam,hell,0,cell_grid_GPT3);
                            a_ngr_w(stnum) = Gn_w;   % [m]
                            a_egr_w(stnum) = Ge_w;   % [m]
                            %a_ngr_w(stnum) = antenna(stnum).gpt3.Gn_w;   % [m]
                            %a_egr_w(stnum) = antenna(stnum).gpt3.Ge_w;   % [m]
                            aprgrd_w = 1/(sin(pi/2-zd)*tan(pi/2-zd)+0.0007)*(a_ngr_w(stnum)*cos(azim)+a_egr_w(stnum)*sin(azim));   % [m]
                    end
                end
                aprgrd = aprgrd_h+aprgrd_w;   % [m]
                a_ngr(stnum) = a_ngr_h(stnum) + a_ngr_w(stnum);   % [m]
                a_egr(stnum) = a_egr_h(stnum) + a_egr_w(stnum);   % [m]

                % Antenna axis offset altitude correction
                aoalt = 0;
                if parameter.vie_mod.aoaltcorr == 1                
                   [psifac] = ao_altcorr(axtyp,ant,LHAe,azim,zd,corz); % Sovers, Fanselow and Jacobs, Eq. 3.202 - 3.206
                   DTSH = 8600; % [m] ~8.6km Dry Troposphere Scale Height
                   aoalt = -zdry*(axkt*c/DTSH)*psifac; %[m]
                end              

                % store in scan
                scan(iSc).stat(stnum).mfw   = mfw;   % used in vie_lsm as partial
                scan(iSc).stat(stnum).zdry  = zdry;
                scan(iSc).stat(stnum).zwet  = zwet;
                scan(iSc).stat(stnum).aoalt   = aoalt; % antenna axis offset altitude correction [m]                    
                scan(iSc).stat(stnum).trop  = ((zdry+aoalt) * mfh + zwet * mfw + aprgrd)/c; % contains the whole (asymmetric) delay [sec]
            end 

            % THERMAL DEFORMATION (Haas et al.,1998; Skurihina 2000)
            if parameter.vie_mod.therm == 1
                thermal = antenna(stnum).thermal;
                temp = scan(iSc).stat(stnum).temp;  % measured temperature [C] or from gpt3
                if isempty(thermal) == 0
                    therm_d = thermdef (...
                        temp,thermal,azim,zd,corz,de,axtyp,offs,aname);
                else
                    therm_d = 0;
                    flagmess.thermal(stnum) = 1;
                end
            else
                therm_d = 0;
            end

            % GRAVITATIONAL DEFORMATION
            if isfield(parameter.vie_mod, 'gravDef') && ...
                     parameter.vie_mod.gravDef == 1 && ...
                     isfield(antenna(stnum), 'gravdef') && ...
                    ~isempty(antenna(stnum).gravdef)                         
                gravdef_data = antenna(stnum).gravdef;
                gravdef_corr = spline(gravdef_data.ez_delay(:,1), gravdef_data.ez_delay(:,2), 90-rad2deg(zd));  % [ps]

                gravdef_corr = gravdef_corr * 1e-12;  % [ps] --> [sec]

                % Exception for ONSALA60
                % elevation independent temperature term
                % delay_temp = (T-19C)*0.47 mm
                if strcmp(strtrim(antenna(stnum).name), 'ONSALA60')
                    temp = scan(iSc).stat(stnum).temp;
                    delay_temp = (temp - 19) * 0.47 * 1e-3;  % [m]
                    delay_temp = delay_temp / c;  % [m] --> [s]

                    gravdef_corr = gravdef_corr + delay_temp;
                end                    
            else
                gravdef_corr = 0;
            end

            % store in scan
            scan(iSc).stat(stnum).axkt    = axkt;
            scan(iSc).stat(stnum).therm   = therm_d;
            scan(iSc).stat(stnum).az      = azim;
            scan(iSc).stat(stnum).zd      = zd;     % corrected for aberration, not for refraction!
            scan(iSc).stat(stnum).daxkt   = daxkt;  % [-] partial derivative AO
            scan(iSc).stat(stnum).gravdef = gravdef_corr;              
        end
    end

    %(8) add geometric part of tropospheric propagation delay
    atm1  = scan(iSc).stat(idStation1).trop;  %[sec]
    atm2  = scan(iSc).stat(idStation2).trop;  %[sec]

    % For the GEOCENTR atm1 becomes -inf!
    % => zenith delay from second station is taken instead as a preliminary solution
    if strcmp(antenna(idStation1).name,'GEOCENTR')
        atm1 = scan(iSc).stat(idStation2).zdry/c;
        scan(iSc).stat(idStation1).trop = atm1;
    elseif strcmp(antenna(idStation2).name,'GEOCENTR')
        atm2 = scan(iSc).stat(idStation1).zdry/c;
        scan(iSc).stat(idStation2).trop = atm2;
    end

    tpd_g = atm1*(rqu*(v2-v1))/c;  
    tau   = tau + tpd_g;         %(eq. 11)

    %(9) total delay
    c_trop = atm2 - atm1;
    tau    = tau + c_trop; %(eq. 12)

    % further corrections:
    % axis offset correction
    c_axis  = scan(iSc).stat(idStation2).axkt - scan(iSc).stat(idStation1).axkt; %[sec]

    % thermal deformation (Attention: Station1 - Station2)
    c_therm = scan(iSc).stat(idStation1).therm - scan(iSc).stat(idStation2).therm; % [sec]

    % gravitational deformation correction
    c_gravdef = scan(iSc).stat(idStation2).gravdef - scan(iSc).stat(idStation1).gravdef; % [sec]

    % add
    tau = c_axis + c_therm + c_gravdef + tau; % [sec]

    % + EXTERNAL IONOSPERIC DELAY +
    if strcmp(parameter.vie_init.iono, 'ext')
        % use iono delay from external file and add it

        % only do the first time
        if ~exist('iondata', 'var')
            if iSc==1
                fprintf('Start loading external ionospheric file\n');
            end
            [iondata, ionFileFoundLog] = load_ionfile(parameter,session);
        end

        % if iono file was found -> apply correction
        if ionFileFoundLog == 1
            % find two stationnames of current observations
            statNameI1      = antenna(scan(iSc).obs(iobs).i1).name;
            statNameI2      = antenna(scan(iSc).obs(iobs).i2).name;
            [dion1, dion2]  = get_iondel(iondata,scan(iSc).tim,statNameI1,statNameI2);
            ionoDel         = dion2 - dion1; %[sec]

            % add to observed delay
            scan(iSc).stat(idStation1).iono  = dion1;
            scan(iSc).stat(idStation2).iono  = dion2;
            scan(iSc).obs(iobs).ionDelext   = ionoDel;
            scan(iSc).obs(iobs).obs         = scan(iSc).obs(iobs).obs - ionoDel;
        else
            warning('Ion. correction not found in external file! Correction set to zero.\n');
            scan(iSc).stat(idStation1).iono  = 0;
            scan(iSc).stat(idStation2).iono  = 0;
            scan(iSc).obs(iobs).ionDelext   = 0;
        end
    end
    %  - EXTERNAL IONOSPERIC DELAY -
end

