% ************************************************************************
%   Description:
%   function to rearrange the data from scan based to source based
%
%   Reference:
%
%   Input:
%       'x'                 (num. of estimates,1)                       all the estimated parameters
%       'first_solution'    (num. of estimates from first solution,1)   num. of estimates from first solution
%       'mi'                (num. of estimates,1)                       formal errors of estimates
%       'na'                (1,1)                                       number of antennas
%       'sum_dj'            (1,number of models)                        total number of estimates for each included model
%       'n_'                structure array                             number of estimates (pwlo or one offset)
%       'mjd0'              (1,1)                                       mjd of the 0:00 UTC of the first day of the session
%       'mjd1'              (1,1)                                       epoch of the first scan in mjd
%       't'                 structure array                             estimation intervals of clk, zwd, ngr, egr, xyz
%       'T'                 structure array                             estimation intervals for dUT1, xpol, ypol, nutdx, and nutdy
%       'opt'               structure array                             (for info. /DOC/opt.doc)
%       'antenna'           structure array                             (for info. /DOC/antenna.doc)
%       'ns_q'              (1,1)                                       number of sources (quasars)
%       'nso'               structure array                             the number of source coordinate offsets for each source
%       'tso'               structure array                             estimation intervals of source coor. for each source
%       'ess'               1/0                                         estimation of parameters from single session
%       'ns_s'              (1,1)                                       number of sources (satellites)
%
%   Output:
%       'x_'                structure array                            (for info. /DOC/x_.doc)
%
%   External calls:
%
%   Coded for VieVS:
%   15 Aug 2009 by Kamil Teke
%
%   Revision:
%   06 Dec 2009 by Kamil Teke: header added
%   13 Oct 2010 by Hana Spicakova: initialisation of x_.ngr.col and x_.egr.col
%   11 Nov 2010 by Hana Spicakova: initialisation of x_.pwclk.col, x_.rqclk.col and x_.zwd.col
%   11 Nov 2010 by Hana Spicakova: indicator for single-session solution added
%   19 May 2011 by Hana Spicakova: initialisation of columns for EOP and st. coordinates, x_.xpol.col, ... and x_.coorx.col,...
%   16 Sep 2011 by Hana Spicakova: initialisation of x_.ngr.mjd and x_.egr.mjd and x_.zwd.mjd
%   19 Jun 2014 by Hana Krasna: sources can be estimated in vie_lsm with NNR condition
%   05 May 2017 by A. Hellerschmied: - general revision
%                                    - Satellite coordinate estimates considered
%   28 Oct 2020 by Johannes Boehm: the values of the first clock solution
%   were wrong in x_, if clock breaks were in the session; this was
%   corrected
%   2025-02-02 by H.Wolf: added Keplerian Element estimation
%   2026-05-27 by H.Wolf: removed parameter ess - one can get it from opt.est_singleses
% ************************************************************************
function [x_] = splitx(x, first_solution, mi, na, sum_dj, n_, mjd0, mjd1, t, T, opt, antenna, ns_q, nso, tso, ns_s, nsat, tsat, ebsl,agnc)

    c = 299792458; % velocity in m/s
    % -------------------------------------------------------------------------
    % DIVIDING THE VECTOR X
    % -------------------------------------------------------------------------
    % clock pw offsets (Split x vector)
    % Units of the estimated clock parameters from first LS solution
    x_.units.pwclk              = 'offset(s) of the clock poly.(s) from the first LS and CPWLO from main LS';
    x_.units.pwclk_first_val    = 'cm';
    x_.units.pwclk_first_mx     = 'cm';
    x_.units.pwclk_first_mjd    = 'epochs of the first scan, break if any, last scan of the antenna';
    x_.units.rclk               = 'rate(s) of the clock poly.(s) from the first LS';
    x_.units.rclk_first_val     = 'cm/day';
    x_.units.rclk_first_mx      = 'cm/day';
    x_.units.rclk_first_mjd     = 'epochs of the first scan, break if any, last scan of the antenna';
    x_.units.qclk               = 'quadratic term(s) of the clock poly.(s) from the first LS';
    x_.units.qclk_first_val     = 'cm/day^2';
    x_.units.qclk_first_mx      = 'cm/day^2';
    x_.units.qclk_first_mjd     = 'epochs of the first scan, break if any, last scan of the antenna';
    
    if opt.first == 1 
        counti = 0; countj = na-1; countk = 2*(na-1);
        
        for istat = 1 : na
            if istat ~= first_solution.ref_st
                if opt.firstclock == 0 || opt.firstclock == 1 || opt.firstclock == 2
                    x_.pwclk(istat).first_mjd = first_solution.breaks(istat).clk; % First clock functions offsets' epochs [mjd]
                    for k = 1 : length(first_solution.breaks(istat).clk) - 1
                        counti = counti + 1;
                        x_.pwclk(istat).first_val(k,:) = first_solution.clk_val(counti,:); % First clock functions offsets [cm]
                        x_.pwclk(istat).first_mx(k,:) = first_solution.clk_mx(counti,:); % formal errors of clock offsets [cm]
                    end
                    countj = countj + length(first_solution.breaks(istat).clk) - 2;
                    countk = countk + length(first_solution.breaks(istat).clk) - 2;
                end
            end
        end
                
        for istat = 1 : na
            if istat ~= first_solution.ref_st
                 if opt.firstclock == 1 || opt.firstclock == 2
                    x_.rclk(istat).first_mjd = first_solution.breaks(istat).clk; % First clock functions offsets' epochs [mjd]
                    for k = 1 : length(first_solution.breaks(istat).clk) - 1
                        countj = countj + 1;
                        x_.rclk(istat).first_val(k,:) = first_solution.clk_val(countj,:); % First clock functions rates [cm/day]
                        x_.rclk(istat).first_mx(k,:) = first_solution.clk_mx(countj,:); % formal errors of clock rates
                    end
                    countk = countk + length(first_solution.breaks(istat).clk) - 2;
                end
            end
        end
                
        for istat = 1 : na
            if istat ~= first_solution.ref_st
                if opt.firstclock == 2
                    x_.qclk(istat).first_mjd = first_solution.breaks(istat).clk; % First clock functions offsets' epochs [mjd]
                    for k = 1 : length(first_solution.breaks(istat).clk) - 1
                        countk = countk + 1;
                        x_.qclk(istat).first_val(k,:) = first_solution.clk_val(countk,:); % First clock functions quadratic terms [cm/day^2]
                        x_.qclk(istat).first_mx(k,:) = first_solution.clk_mx(countk,:); % formal errors of clock quadratic terms
                    end
                end
            end
        end
        
    end 
    
    % Units of the estimated CPWLO clock parameters from main LS solution
    x_.units.pwclk_val = 'from main LS: CPWLO of clocks in cm';
    x_.units.pwclk_mx = 'formal errors of CPWLO of clocks in cm';
    x_.units.pwclk_mjd = 'epochs of CPWLO of clocks in mjd';
    x_.units.pwclk_col = 'column and row numbers of CPWLO clock estimates in N and b';
    
    sumclk = 0;
    for istat = 1 : na
        x_.pwclk(istat).col=[]; % 11Nov10 hana
        for ioffset = 1 : n_(istat).clk
            if opt.est_singleses
                x_.pwclk(istat).val(ioffset,:) = x(sum_dj(1) + sumclk + ioffset,:); % [cm] - estimated VALue
                x_.pwclk(istat).mx(ioffset,:) = mi(sum_dj(1) + sumclk + ioffset,:); % [cm] - STD deviation of the estimate
            end
            x_.pwclk(istat).mjd(ioffset) = t(istat).clk(ioffset)/(24*60) + mjd0; % [MJD] - time of the estimate
            x_.pwclk(istat).col(ioffset) = sum_dj(1) + sumclk + ioffset; % [1] - COLumn of the estimate in A or N
        end
        sumclk = sumclk + n_(istat).clk;
    end
    
    % -------------------------------------------------------------------------
    % clock rate & quadratic terms (Split x vector)
    % Units of the estimated clock parameters from main LS solution
    x_.units.rqclk = 'from main LS: rate and quadratic terms of the clock poly.';
    x_.units.rqclk_val = '(1)cm/day and (2)cm/day^2';
    x_.units.rqclk_mx = '(1)cm/day and (2)cm/day^2';
    x_.units.rqclk_col = 'column and row numbers of the clock rate and quadratic terms in N and b';
    
    sumqclk = 0;
    for istat = 1 : na
        x_.rqclk(istat).col=[]; % 11Nov10 hana
        for irate = 1 : n_(istat).qclk
            if opt.est_singleses
                x_.rqclk(istat).val(irate,:) = x(sum_dj(2) + sumqclk + irate,:); % estimated VALue
                x_.rqclk(istat).mx(irate,:) = mi(sum_dj(2) + sumqclk + irate,:); % STD deviation of the estimate
            end
            x_.rqclk(istat).col(irate) = sum_dj(2) + sumqclk + irate; % [1] - COLumn of the estimate in A or N
        end
        sumqclk = sumqclk + n_(istat).qclk;
    end
    
    % -------------------------------------------------------------------------
    % zwd pw offsets (Split x vector)
    % Units of the estimated ZWD
    x_.units.zwd_val = 'ZWD CPWLO in cm';
    x_.units.zwd_mx = 'formal errors of ZWD CPWLO in cm';
    x_.units.zwd_mjd = 'epochs of ZWD CPWLO in mjd';
    x_.units.zwd_col = 'column and row numbers of ZWD in N and b';
    
    sumzwd = 0;
    for istat = 1:na
        x_.zwd(istat).col=[]; % 11Nov10 hana
        x_.zwd(istat).mjd=[]; % 16Sep11 hana
        x_.antenna(istat).name = antenna(istat).name;
        for ioffset = 1 : n_(istat).zwd
            if opt.est_singleses
                x_.zwd(istat).val(ioffset,:) = x(sum_dj(3) + sumzwd + ioffset,:); % [cm] estimated VALue
                x_.zwd(istat).mx(ioffset,:) = mi(sum_dj(3) + sumzwd + ioffset,:); % [cm] STD deviation of the estimate
            end
            x_.zwd(istat).mjd(ioffset) = t(istat).zwd(ioffset)/(24*60) + mjd0; % [MJD] - time of the estimate
            x_.zwd(istat).col(ioffset) = sum_dj(3) + sumzwd + ioffset; % [1] - COLumn of the estimate in A or N
        end
        sumzwd = sumzwd + n_(istat).zwd;
    end
    
    % -------------------------------------------------------------------------
    % north gradients pw offsets (Split x vector)
    % Units of the estimated troposphere north gradients
    x_.units.ngr_val = 'CPWLO of troposphere north gradients in cm';
    
    sumngr = 0;
    for istat = 1:na
        x_.ngr(istat).col=[]; % 13Oct10 hana
        x_.ngr(istat).mjd=[]; % 16Sep11 hana
        for ioffset = 1 : n_(istat).ngr
            if opt.est_singleses
                x_.ngr(istat).val(ioffset,:) = x(sum_dj(4) + sumngr + ioffset,:); % [cm] estimated VALue
                x_.ngr(istat).mx(ioffset,:) = mi(sum_dj(4) + sumngr + ioffset,:);  % [cm] STD deviation of the estimate
            end
            x_.ngr(istat).mjd(ioffset) = t(istat).ngr(ioffset)/(24*60) + mjd0; % [MJD] - time of the estimate
            x_.ngr(istat).col(ioffset) = sum_dj(4) + sumngr + ioffset; % [1] - COLumn of the estimate in A or N
        end
        sumngr = sumngr + n_(istat).ngr;
    end
    
    % -------------------------------------------------------------------------
    % east gradients pw offsets (Split x vector)
    % Units of the estimated troposphere east gradients
    x_.units.egr_val = 'CPWLO of troposphere east gradients in cm';
    
    sumegr = 0;
    for istat = 1:na
        x_.egr(istat).col=[]; % 13Oct10 hana
        x_.egr(istat).mjd=[]; % 16Sep11 hana
        for ioffset = 1 : n_(istat).egr
            if opt.est_singleses
                x_.egr(istat).val(ioffset,:) = x(sum_dj(5) + sumegr + ioffset,:); % [cm] estimated VALue
                x_.egr(istat).mx(ioffset,:) = mi(sum_dj(5) + sumegr + ioffset,:);  % [cm] STD deviation of the estimate
            end
            x_.egr(istat).mjd(ioffset) = t(istat).egr(ioffset)/(24*60) + mjd0; % [MJD] - time of the estimate
            x_.egr(istat).col(ioffset) = sum_dj(5) + sumegr + ioffset; % [1] - COLumn of the estimate in A or N
        end
        sumegr = sumegr + n_(istat).egr;
    end
    
    % -------------------------------------------------------------------------
    % XPOL pw offsets (Split x vector)
    % Units of xpol
    x_.units.xpol_val = 'mas';
    
    x_.xpol.col=[]; %19May11 hana
    x_.xpol.mjd=[]; %19May11 hana
    
    for ioffset = 1 : size(T.xpol,2)
        if opt.est_singleses
            x_.xpol.val(ioffset,:) = x(sum_dj(6) + ioffset,:); % [mas] estimated VALue
            x_.xpol.mx(ioffset,:) = mi(sum_dj(6) + ioffset,:); % [mas] STD deviation of the estimate
        end
        x_.xpol.mjd(ioffset) = T.xpol(ioffset)/(24*60) + mjd0; % [MJD] - time of the estimate
        x_.xpol.col(ioffset) = sum_dj(6) + ioffset; % [1] - COLumn of the estimate in A or N
    end
    
    % -------------------------------------------------------------------------
    % YPOL pw offsets (Split x vector)
    % Units of ypol
    x_.units.ypol_val = 'mas';
    
    x_.ypol.col=[]; %19May11 hana
    x_.ypol.mjd=[]; %19May11 hana
    
    for ioffset = 1 : size(T.ypol,2)
        if opt.est_singleses
            x_.ypol.val(ioffset,:) = x(sum_dj(7) + ioffset,:); % [mas] estimated VALue
            x_.ypol.mx(ioffset,:) = mi(sum_dj(7) + ioffset,:); % [mas] STD deviation of the estimate
        end
        x_.ypol.mjd(ioffset) = T.ypol(ioffset)/(24*60) + mjd0; % [MJD] - time of the estimate
        x_.ypol.col(ioffset) = sum_dj(7) + ioffset; % [1] - COLumn of the estimate in A or N
    end
    
    % -------------------------------------------------------------------------
    % dUT1 pw offsets (Split x vector)
    % Units of dut1
    x_.units.dut1_val = 'ms';
    x_.units.dut1_mx = 'ms';
    
    x_.dut1.col=[]; %19May11 hana
    x_.dut1.mjd=[]; %19May11 hana
    
    for ioffset = 1 : size(T.dut1,2)
        if opt.est_singleses
            x_.dut1.val(ioffset,:) = x(sum_dj(8) + ioffset,:)/15; % [ms] estimated VALue
            x_.dut1.mx(ioffset,:) = mi(sum_dj(8) + ioffset,:)/15; % [ms] STD deviation of the estimate
        end
        x_.dut1.mjd(ioffset) = T.dut1(ioffset)/(24*60) + mjd0; % [MJD] - time of the estimate
        x_.dut1.col(ioffset) = sum_dj(8) + ioffset; % [1] - COLumn of the estimate in A or N
    end
    
    % -------------------------------------------------------------------------
    % Celestial Pole coordinates dx (DEPS) pw offsets (Split x vector)
    % Units of nutation offsets
    x_.units.nutdx_val = 'mas';
    
    x_.nutdx.col=[]; %19May11 hana
    x_.nutdx.mjd=[]; %19May11 hana
    
    for ioffset = 1 : size(T.nutdx,2)
        if opt.est_singleses
            x_.nutdx.val(ioffset,:) = x(sum_dj(9) + ioffset,:); % [mas] estimated VALue
            x_.nutdx.mx(ioffset,:) = mi(sum_dj(9) + ioffset,:); % [mas] STD deviation of the estimate
        end
        x_.nutdx.mjd(ioffset) = T.nutdx(ioffset)/(24*60) + mjd0; % [MJD] - time of the estimate
        x_.nutdx.col(ioffset) = sum_dj(9) + ioffset; % [1] - COLumn of the estimate in A or N
    end
    
    % -------------------------------------------------------------------------
    % Celestial Pole coordinates dx (PSI) pw offsets (Split x vector)
    % Units of nutation offsets
    x_.units.nutdy_val = 'mas';
    
    x_.nutdy.col=[]; %19May11 hana
    x_.nutdy.mjd=[]; %19May11 hana
    
    for ioffset = 1 : size(T.nutdy,2)
        if opt.est_singleses
            x_.nutdy.val(ioffset,:) = x(sum_dj(10) + ioffset,:); % [mas] estimated VALue
            x_.nutdy.mx(ioffset,:) = mi(sum_dj(10) + ioffset,:); % [mas] STD deviation of the estimate
        end
        x_.nutdy.mjd(ioffset) = T.nutdy(ioffset)/(24*60) + mjd0; % [MJD] - time of the estimate
        x_.nutdy.col(ioffset) = sum_dj(10) + ioffset; % [1] - COLumn of the estimate in A or N
    end
    
    % Source coordinates
    % Units of CPWLO source coordinates
    if opt.pw_sou == 1
        x_.units.soura_val = 'source CRF CPWLO right ascension coor. in ms';
        x_.units.soura_mx = 'ms';
        x_.units.soude_val = 'source CRF CPWLO declination coor. in mas';
        x_.units.soude_mx = 'mas';
    end
    if opt.est_sourceNNR==1
        x_.units.soura_val = 'source CRF RA estimate (with NNR condition) in ms';
        x_.units.soura_mx = 'ms';
        x_.units.soude_val = 'source CRF De estimate (with NNR condition) in mas';
        x_.units.soude_mx = 'mas';
    end
    
    
    
    for jsou = 1 : ns_q
        x_.soura(jsou).val = [];
        x_.soura(jsou).mjd = [];
        x_.soura(jsou).mx = [];
        x_.soura(jsou).col = [];
        x_.soura(jsou).inNNR = [];
        x_.soura(jsou).crade = [];
        x_.soude(jsou).val = [];
        x_.soude(jsou).mjd = [];
        x_.soude(jsou).mx = [];
        x_.soude(jsou).col = [];
    end
    
    
    if opt.pw_sou == 1
        sumsou = 0; jsou = 0;
        for isou = 1 : ns_q
            if opt.source(isou).rade_inc == 1
                jsou = jsou + 1;
                % -------------------------------------------------------------------------
                % source coordinates (right ascension) pw offsets (Split x vector)
                x_.soura(jsou).name = opt.source(isou).name;
                if opt.est_singleses
                    for ioffset = 1 : nso(isou).sources
                        x_.soura(jsou).val(ioffset,:) = x(sum_dj(11) + sumsou + ioffset,:)/15; % [ms] estimated value
                        %x_.soura(jsou).mjd(ioffset)  = tso(jsou).sources(ioffset)/(24*60) + mjd0; % [MJD] - time of the estimate; HK: doesn't work after the orbit update
                        x_.soura(jsou).mjd(ioffset)  = tso(isou).sources(ioffset)/(24*60) + mjd0; % [MJD] - time of the estimate; HK: changed
                        x_.soura(jsou).mx(ioffset,:) = mi(sum_dj(11) + sumsou + ioffset,:)/15; % [ms] std. dev. of the estimate
                        x_.soura(jsou).col(ioffset)  = sum_dj(11) + sumsou + ioffset; % [1] - COLumn of the estimate in A or N
                        x_.soura(jsou).crade(ioffset)  =agnc(sumsou + ioffset); % correlation
                    end
                end
                % -------------------------------------------------------------------------
                % source coordinates (declination) pw offsets (Split x vector)
                x_.soude(jsou).name = opt.source(isou).name;
                if opt.est_singleses
                    for ioffset = 1 : nso(isou).sources
                        x_.soude(jsou).val(ioffset,:) = x(sum_dj(12) + sumsou + ioffset,:); % [mas] estimated value
                        %x_.soude(jsou).mjd(ioffset) = tso(jsou).sources(ioffset)/(24*60) + mjd0; % [MJD] - time of the estimate; HK: doesn't work after the orbit update
                        x_.soude(jsou).mjd(ioffset) = tso(isou).sources(ioffset)/(24*60) + mjd0; % [MJD] - time of the estimate; HK: changed
                        x_.soude(jsou).mx(ioffset,:) = mi(sum_dj(12) + sumsou + ioffset,:); % [mas] std. dev. of the estimate;
                        x_.soude(jsou).col(ioffset) = sum_dj(12) + sumsou + ioffset; % [1] - COLumn of the estimate in A or N
                    end
                end
                sumsou = sumsou + nso(isou).sources;
            end
        end
    end
    
    % sources estimated with NNR condition
    if opt.est_singleses
        if opt.est_sourceNNR==1
            for jsou = 1 : ns_q
                x_.soura(jsou).val = x(sum_dj(11) + jsou ,:)/15; % [ms] estimated value
                x_.soura(jsou).mjd = ceil(mjd1); % mjd1 : midnight
                x_.soura(jsou).mx = mi(sum_dj(11) + jsou ,:)/15; % [ms] std. dev. of the estimate
                x_.soura(jsou).col = sum_dj(11) + jsou ; % [1] - COLumn of the estimate in A or N
                x_.soura(jsou).inNNR = opt.source(jsou).nnr_inc; % included in NNR 1/0
                x_.soura(jsou).crade =agnc(jsou); % correlation
                
                x_.soude(jsou).val = x(sum_dj(12) + jsou ,:); % [mas] estimated value
                x_.soude(jsou).mjd = ceil(mjd1); % mjd1 : midnight
                x_.soude(jsou).mx = mi(sum_dj(12) + jsou ,:); % [mas] std. dev. of the estimate
                x_.soude(jsou).col = sum_dj(12) + jsou ; % [1] - COLumn of the estimate in A or N
                x_.soude(jsou).inNNR = opt.source(jsou).nnr_inc; % included in NNR 1/0
            end
        end
    end
    
    % Antenna coordinates
    % Units of antenna coordinates
    x_.units.coorx_val = 'antenna TRF X coor. estimate(s) in cm (using all observations)';
    x_.units.coorx_mx = 'cm';
    x_.units.coory_val = 'antenna TRF Y coor. estimate(s) in cm (using all observations)';
    x_.units.coory_mx = 'cm';
    x_.units.coorz_val = 'antenna TRF Z coor. estimate(s) in cm (using all observations)';
    x_.units.coorz_mx = 'cm';

    x_.units.coorx_sat_val = 'antenna TRF X coor. estimate(s) in cm (using satellite obs. only)';
    x_.units.coorx_sat_mx = 'cm';
    x_.units.coory_sat_val = 'antenna TRF Y coor. estimate(s) in cm (using satellite obs. only)';
    x_.units.coory_sat_mx = 'cm';
    x_.units.coorz_sat_val = 'antenna TRF Z coor. estimate(s) in cm (using satellite obs. only)';
    x_.units.coorz_sat_mx = 'cm';

    x_.units.coorx_qu_val = 'antenna TRF X coor. estimate(s) in cm (using quasar obs. only)';
    x_.units.coorx_qu_mx = 'cm';
    x_.units.coory_qu_val = 'antenna TRF Y coor. estimate(s) in cm (using quasar obs. only)';
    x_.units.coory_qu_mx = 'cm';
    x_.units.coorz_qu_val = 'antenna TRF Z coor. estimate(s) in cm (using quasar obs. only)';
    x_.units.coorz_qu_mx = 'cm';
    
    switch true
        case opt.stc_sat
            base_offs = [27, 28, 29];
            flds      = {'coorx_sat', 'coory_sat', 'coorz_sat'};
            n_cfg     = 1;
        case opt.stc_qu
            base_offs = [30, 31, 32];
            flds      = {'coorx_qu', 'coory_qu', 'coorz_qu'};
            n_cfg     = 1;
        case opt.stc_all
            base_offs = [13, 14, 15];
            flds      = {'coorx', 'coory', 'coorz'};
            n_cfg     = 1;
        case opt.stc_qs
            base_offs = {[27, 30], [28, 31], [29, 32]};
            flds      = {{'coorx_sat', 'coorx_qu'}, {'coory_sat', 'coory_qu'}, {'coorz_sat', 'coorz_qu'}};
            n_cfg     = 2;
    end

    sumxyz = 0;
    for istat = 1:na
        n_xyz = n_(istat).xyz;
        if n_xyz == 0, continue; end
        
        if opt.pw_stc == 0
            mjds = ceil(mjd1) * ones(n_xyz, 1);
        else
            mjds = t(istat).xyz(1:n_xyz) / (24*60) + mjd0;
        end
        
        for c = 1:3
            if n_cfg == 1
                cur_base = base_offs(c);
                cur_names = {flds{c}};
            else % stc_qs
                cur_base = base_offs{c};
                cur_names = flds{c};
            end
            
            for k = 1:length(cur_names)
                fname = cur_names{k};
                idx = sum_dj(cur_base(k)) + sumxyz + (1:n_xyz);
                
                x_.(fname)(istat).col = idx';
                x_.(fname)(istat).mjd = mjds';
                x_.(fname)(istat).val = [];
                x_.(fname)(istat).mx  = [];
                
                if opt.est_singleses
                    x_.(fname)(istat).val = x(idx, :);
                    x_.(fname)(istat).mx  = mi(idx, :);
                end
            end
        end
        sumxyz = sumxyz + n_xyz;
    end
    
    % -------------------------------------------------------------------------
    % Satellite coordinate offsets (PWL) 
    % -------------------------------------------------------------------------
    % Units:
    frame = opt.SatPos.sat_pos_est_ref_frame;
    switch frame
        case {'gcrf', 'trf'}
            compNames = {'X', 'Y', 'Z'};
        case 'rsw'
            compNames = {'R', 'S', 'W'};
        case 'ntw'
            compNames = {'N', 'T', 'W'};
    end

    ns_s_tmp = max(ns_s, 1);
    compBaseOff = [16, 17, 18];
    for iComp = 1:3
        compStr = compNames{iComp};
        fName   = sprintf('sat_pos%s', compStr); 
        
        x_.units.(fName).val = sprintf('satellite coor. (%s) estimate(s) in cm, ref.-frame: %s', ...
                                       compStr, frame);
        x_.units.(fName).mx  = 'cm';  
        x_.(fName)(ns_s_tmp).name = [];
        x_.(fName)(ns_s_tmp).val = [];
        x_.(fName)(ns_s_tmp).mjd = [];
        x_.(fName)(ns_s_tmp).mx  = [];
        x_.(fName)(ns_s_tmp).col = [];
        compFields{iComp} = fName;
    end
      
    if opt.SatPos.pw_sat == 1
        sum_sat = 0;  % Sum of satelite pos estimates
        isat2 = 0;  % sat index
        
        % Loop over all satellites:
        for isat = 1 : ns_s
            % If the pos. of the current satellite was estimated
            if opt.satellite(isat).SatPos.pos_inc == 1
                isat2 = isat2 + 1;
                nPos  = nsat(isat2).pos;
                if nPos == 0, continue; end

                mjds = tsat(isat).pos(1:nPos)/(24*60) + mjd0;
                for c = 1:3
                    fName = compFields{c};
                    off   = sum_dj(compBaseOff(c)) + sum_sat;    
                    idx = off + (1:nPos);
                    
                    x_.(fName)(isat2).name = opt.satellite(isat).name;
                    x_.(fName)(isat2).val  = x(idx, :);
                    x_.(fName)(isat2).mjd  = mjds';
                    x_.(fName)(isat2).mx   = mi(idx, :);
                    x_.(fName)(isat2).col  = idx';
                end                
                sum_sat = sum_sat + nsat(isat2).pos;
             end 
         end 
    end
    
    x_.units.scale = 'correction to the scale [ppb]';
    x_.scale.col = []; %[-]
    x_.scale.val = []; %[-]
    x_.scale.mx = []; %[-]
    x_.scale.mjd = []; %[-]
    
    x_.units.bdclko = 'baseline-dependent clock offset [cm]';
    x_.bdclko.val = []; % cm
    x_.bdclko.mx = [];
    x_.bdclko.mjd = [];
    x_.bdclko.col = [];
    x_.bdclko.st1 = [];
    x_.bdclko.st2 = [];
    x_.bdclko.namest1 = [''];
    x_.bdclko.namest2 = [''];
    
    % Correction to the scale factor
    if opt.est_scale==1
        if opt.est_singleses
            x_.scale.val = x(sum_dj(20)) /c/100 *1e9; %[ppb]
            x_.scale.mx = mi(sum_dj(20)) /c/100 *1e9; %[ppb]
            x_.scale.col = sum_dj(20);
            x_.scale.mjd = ceil(mjd1); % mjd1 : midnight
        else
            x_.scale.col = sum_dj(20);
            x_.scale.mjd = ceil(mjd1); % mjd1 : midnight
        end
    end
    
    % Baseline-dependent clock offset
    if opt.est_bdco==1
        nbas = sum_dj(21) - sum_dj(20);
        if opt.est_singleses
            for i=1:nbas
                x_.bdclko(i).val = x(sum_dj(20)+i); % cm
                x_.bdclko(i).mx = mi(sum_dj(20)+i);
                x_.bdclko(i).col = sum_dj(20)+i;
                x_.bdclko(i).st1 = ebsl(i,1);
                x_.bdclko(i).st2 = ebsl(i,2);
                x_.bdclko(i).namest1 = antenna(ebsl(i,1)).name;
                x_.bdclko(i).namest2 = antenna(ebsl(i,2)).name;
                x_.bdclko(i).mjd = ceil(mjd1); 
            end
        else
            for i=1:nbas
                x_.bdclko(i).col = sum_dj(20)+i;
            end
        end
    end
    
    
    % -------------------------------------------------------------------------
    % Keplerian Elements (PWL)
    % -------------------------------------------------------------------------
    % Units
    x_.units.ORB.sma_val   = 'orbital element semimajor axis a estimate(s) in cm'; 
    x_.units.ORB.sma_mx    = 'cm';
    x_.units.ORB.ecc_val   = 'orbital element eccentricity estimate(s) unitless'; 
    x_.units.ORB.ecc_mx    = '[]'; 
    x_.units.ORB.inc_val   = 'orbital element inclination estimate(s) in mas'; 
    x_.units.ORB.inc_mx    = 'mas';  
    x_.units.ORB.raan_val   = 'orbital element RAAN estimate(s) in mas'; 
    x_.units.ORB.raan_mx    = 'mas';
    x_.units.ORB.argp_val   = 'orbital element argument of periapsis estimate(s) in mas'; 
    x_.units.ORB.argp_mx    = 'mas'; 
    x_.units.ORB.argl_val   = 'orbital element argument of latitude estimate(s) in mas'; 
    x_.units.ORB.argl_mx    = 'mas'; 
            
    iMat = 20;
    orbNames = ["sma"; "ecc"; "inc"; "raan"; "argp"; "argl"];

    for iorb = 1:length(opt.ORB.params)
        x_.('ORB').(orbNames(iorb)).name = [];
        x_.('ORB').(orbNames(iorb))(ns_s_tmp).val   = [];
        x_.('ORB').(orbNames(iorb))(ns_s_tmp).mjd   = [];
        x_.('ORB').(orbNames(iorb))(ns_s_tmp).mx    = [];
        x_.('ORB').(orbNames(iorb))(ns_s_tmp).col   = [];

        if opt.ORB.estORB && opt.ORB.params(iorb).estimate
            sum_sat = 0;  % Sum of satellite position estimates
            for i_sat = 1:ns_s
                x_.('ORB').(orbNames(iorb))(i_sat).name = opt.satellite(i_sat).name;
                nsat_count = nsat(i_sat).('orb' + string(iorb));
                offsets = sum_dj(iMat + iorb) + sum_sat + (1:nsat_count);
    
                x_.('ORB').(orbNames(iorb))(i_sat).val(1:nsat_count, :) = x(offsets, :);
                x_.('ORB').(orbNames(iorb))(i_sat).mjd(1:nsat_count) = tsat(1).('orb' + string(iorb))(1:nsat_count) / (24 * 60) + mjd0;
                x_.('ORB').(orbNames(iorb))(i_sat).mx(1:nsat_count, :) = mi(offsets, :);
                x_.('ORB').(orbNames(iorb))(i_sat).col(1:nsat_count) = offsets';
    
                sum_sat = sum_sat + nsat_count;
            end
        end
    end

    % -------------------------------------------------------------------------
    % Dynamical Parameters (PWL)
    % -------------------------------------------------------------------------
    % Units
    x_.units.SRP_val   = 'direct rpr estimate(s)'; 
    x_.units.SRP_mx    = 'nm/s^2';
        
    iMat = 32;
    srpNames = ["D0"; "Y0"; "B0"; "DC"; "YC"; "BC"; "DS"; "YS"; "BS"];
      
    for isrp = 1:length(opt.SRP.params)
        x_.('SRP').(srpNames(isrp))(ns_s_tmp).name = [];
        x_.('SRP').(srpNames(isrp))(ns_s_tmp).val   = [];
        x_.('SRP').(srpNames(isrp))(ns_s_tmp).mjd   = [];
        x_.('SRP').(srpNames(isrp))(ns_s_tmp).mx    = [];
        x_.('SRP').(srpNames(isrp))(ns_s_tmp).col   = [];

        if opt.SRP.estSRP && opt.SRP.params(isrp).estimate 
            sum_sat = 0;
            for isat = 1:ns_s
                x_.('SRP').(srpNames(isrp))(isat).name = opt.satellite(isat).name;
                nsat_count = nsat.('srp' + string(isrp));
                offsets = sum_dj(iMat + isrp) + sum_sat + (1:nsat_count);

                x_.('SRP').(srpNames(isrp))(isat).val(1:nsat_count, :) = x(offsets, :);
                x_.('SRP').(srpNames(isrp))(isat).mjd(1:nsat_count) = tsat(1).('srp' + string(isrp))(1:nsat_count) / (24 * 60) + mjd0;
                x_.('SRP').(srpNames(isrp))(isat).mx(1:nsat_count, :) = mi(offsets, :);
                x_.('SRP').(srpNames(isrp))(isat).col(1:nsat_count) = offsets';
                sum_sat = sum_sat + nsat_count;
            end
        end
    end  
end