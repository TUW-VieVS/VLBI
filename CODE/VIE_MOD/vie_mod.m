% ************************************************************************
%   Description:
%   In Vie_mod the theoretical time delay for every baseline is computed,
%   as well as its partial derivatives. Therefore, the source and station
%   coordinates are treated, the ephemerides, the Earth Orientation and
%   various models for station corrections. The program is arranged in a
%   common part (for all observation epochs) and loops over the scans,
%   antennas and baselines.
%
%   References:
%    - IERS Conventions 2000 plus updates (as of 5 Nov 2009)
%    - Sovers et al. 1998 "Astrometry and geodesy with radio
%       interferometry: experiments, models, results"
%    - see subfunctions for details
%
%   Input:
%      ---
%      the program is called from vievs.m
%      Data input is done using the structure arrays:
%               - currentf_so, currentf_an, currentf_sc, currentf_pa
%       from the WORK directory
%
%   Output:
%      the output data is stored in the structure arrays
%           - currentf_so, currentf_an, currentf_sc, currentf_pa in the
%               WORK directory and
%           - name of session_sources, _antenna, _scan, _parameter in
%               .../DATA/LEVEL1/
%
%   External calls:
%   constants.m, tai_utc.m, load_eop.m, rad2as.m, tver2000.m, lagint4v.m,
%   eophf.m, load_eph.m, trs2crs.m, sourcevev.m, corr_ant.m, matthew.m,
%   ctocean.m, ctatmos12.m, ctatmo.m, ctpole.m, xyz2ell.m, locsource.m,
%   thermdef.m, axis_stat.m, parlovshi.m, grav_delay.m,
%   call_spline_4.m, calc_met_data
%
%   Coded for VieVS:
%   05 Nov 2009 by Lucia Plank
%
%   Revision:
%   17 Dec 2009 by Lucia Plank: replace xyz2ellip.m with xyz2ell.m
%   24 Feb 2010 by Lucia Plank and Tobias Nilsson: file paths for input
%      and output are modified user defined sub-directory for input and
%      output is available.
%   01 Mar 2010 by Lucia Plank: changes in scan array --> ephemerides are
%      not included any more. Now the ephemerides are stored in
%      ../EPHEM/ephem_sessionname.mat . If this file already exists, the
%      ephemerides are not calculated again.
%   05 Mar 2010 by Lucia Plank: Ephemerides JPL_421 implemented.
%   12 Mar 2010 by Lucia Plank: Antenna eccentricities added.
%   30 Mar 2010 by Lucia Plank: scan.tim
%   25 May 2010 by Hana Spicakova: new partials for Love and Shida numbers
%   26 May 2010 by Lucia Plank: Aberration added;
%   06 Aug 2010 by Lucia Plank: grav_delay new
%   13 Oct 2010 by Hana Spicakova: model for conventional mean pole added
%      according to IERS Conv. 2010 (cubic model over the period 1976-2010)
%   12 Nov 2010 by Hana Spicakova: big update concerning frequency
%       dependent Love/Shida numbers (see matthew.m)
%   12 Nov 2010 by Hana Spicakova: partial derivatives wrt FCN period are
%       added (IERS03 Ch6 Eq.6) (see matthew.m)
%   08 Feb 2011 by Hana Spicakova: models for a priori troposp. gradient
%       models added (apg.m and dao.m)
%   08 Apr 2011 by Matthias Madzak: added possibility of using external
%       tropospheric files
%   18 Apr 2011 by Tobias Nilsson: Changes input/output variables
%   04 May 2011 by Matthias Madzak: Changed loop over baselines to add
%       ionospheric delay from TEC map when button in vie_mod-Gui ("use
%       external iono") is ticked.
%   05 May 2011 by Matthias Madzak: Changed welcome message to be
%       consistent with vie_init and vie_lsm.
%   17 May 2011 by Hana Spicakova: a priori tr. gradients are stored into
%       antenna structure
%   19 Jul 2011 by Hana Spicakova: displacement of the points due to tidal
%       ocean loading is computed according to IERS Conv. 2010, Ch 7.1.2. The
%       subroutine HARDISP.F provided by D. Agnew was translated into Matlab
%       and rearranged
%   18 Aug 2011 by Lucia Plank: major revision, harmonization of recent
%       changes
%   10 Nov 2011 by Hana Spicakova: popup menu for APL series
%   10 Nov 2011 by Hana Spicakova: ocean pole tide loading added (IERS2010)
%   27 Mar 2012 by Hana Spicakova: hydrology loading
%   21 May 2012 by Lucia Plank: major revision, harmonization of recent
%       changes; removing the calculation of external trop. & ion. delays from
%       vie_mod to subroutines.
%   19 Jun 2012 by Hana Krasna: 
%          partial derivatives of Love and Shida numbers, FCN period from
%          solid Earth tides, partials w.r.t. acceleration of SSB, partials
%          w.r.t. amplitudes of seasonal variation in station positions
%   25 Sep 2012 by Hana Krasna: partials w.r.t. relativistic parameter
%          gamma added  (only Sun's contribution)
%   02 Nov 2012 by Lucia Plank: enable tidal UT for linear interpolation
%   08 Nov 2012 by Hana Krasna: changed to superstation file (ocean tidal loading,
%                  atmosphere tidal loading, ocean pole tide loading,
%                  thermal deformation coeeficients, axis offset)
%   08 Nov 2012 by Hana Krasna: added GPT2 mapping function
%   04 Oct 2013 by Hana Krasna: partials for antenna axis offset added
%   15 Oct 2013 by Hana Krasna: parameters for internal version moved into
%          vie_setup
%   05 Dec 2013 by Hana Krasna: APL regression coefficients added
%         
%   23 Jan 2014 by Hana Krasna: GIA uplift rates as apriori added 
%   05 Feb 2014 by Lucia Plank: source structure correction + Jet files
%   12 Sep 2014 by Lucia Plank: source structure correction: IVS/IERS
%           source name check
%   16 Dec 2014 by Daniel Landskron: nGrad and eGrad deleted (not used)
%   21 Aug 2015 by Sigrid Boehm: implementation of partial set-up for
%           high-frequency tidal ERP variations.
%   01 Sep 2015 by Daniel Landskron: get_trpdel gets a further input
%          parameter "sourceNames"
%   21 Oct 2015 by Hana Krasna: partials w.r.t. amplitudes of seasonal variation
%          in station positions, APL regression coefficients, pole tide Love numbers,
%          GIA uplift rates as apriori added.
%   03 Nov 2015 by Armin Hofmeister: In case trp data should be used, but it is missing then
%          display message about type of mapping function used instead of trp data.
%          Use logical variable trp_used instead of multiple equal strcmp calls.
%   22 Dec 2015 by Armin Hofmeister: changes due to the removal of the output parameter
%          "vahabsRaytracingFiles" in function load_trpfile.m
%   26 Apr 2016 by Daniel Landskron: lagint4v/lagint4 changed to spline interpolation (except for EOP's)
% ************************************************************************
%   29 Apr 2016 by M. Madzak: Added post-seismic deformation, e.g. for ITRF2014
%   11 Jun 2016 by Hana Krasna: implementation of partials for the
%   FCN period and amplitudes
%   23 Jun 2016 by A. Girdiuk: bug-fix in loadings
%   26 Jun 2016 by A. Hellerschmied: Added call of function "calc_met_data" to calculate meterologic data with the GPT and GPT2 models. This was done befor in vie_init.
%   06 Jul 2016 by A. Hellerschmied: Notifications are written to the CW if met. data from backup models was used.
%   17 Jul 2016 by A. Hellerschmied: Content of "sources" structure (natural sources, quasars) is now stored in the sub-structure "sources.q"
%   08 Aug 2016 by A. Girdiuk: new ephemeris (jpl_430) and hydrology loading (ERAHYD) provider are added
%   31 Aug 2016 by A. Girdiuk: Minor change: masses of celestial bodies are stored in the ephem structure also for jpl_430 data.
%   31 Aug 2016 by A. Hellerschmied: - Support of near field targets:
%                                       - New source type supported: Spacecraft
%                                       - Near field model added (Iterative solution of the light time equation, taken from vie_mod_tie.m by L. Plank)
%                                    - Bug fix: There was an error at loading existing ephemerids files
%   02 Sep 2016 by A. Hellerschmied: Bug fixes and minor changes subsequent to the update from 31 Aug 2016
%                                    - NOTE: "rq" and "rqu" are probably not correctly calculated for spacecrafts!!!!!
%   02 Sep 2016 by A. Hellerschmied: Changes to be able to calculate "geocentric" a delays:
%                                     - first station of baseline has to be "GEOCENTR"; "GEOCENTR" needs to be defined in the superstaion file
%                                     - For "GEOCENTR" some corrections (e.g. station displacement) are overridden 
%   06 Oct 2016 by Armin Hofmeister: adapt output message of backup trpdata in case no trp file is found
%   14 Oct 2016 by A. Hellerschmied: Revised ion. delay correction from external files
%   22 Nov 2016 by A. Hellerschmied: Added posibility to use SC velocities if provided in the sources structure.
%   05 Dec 2016 by A. Hellerschmied: - Instead of MJD epochs "seconds since reference epoch" (given by integer MJD and integer seconds) is now used to interpolate S.C. pos. and vel. for a priori delays.
%                                      This was required, because MJD provided limited time resolution only (~1 musec signifiancy limit)
%                                    - Number of iteration of S.C. vel. and pos. is adaptive now (set: "ddt_threshold" and "max_iterations")
%   18 Jan 2017 by D. Landskron: a priori gradient section extended
%   20 Jan 2017 by D. Landskron: also ray-tracing .trp-files as input enabled
%   23 Jan 2017 by D. Landskron: GPT2 changed to GPT3, GPT removed
%   07 Jan 2017 by A. Hellerschmied: Several bug-fixes in subsequent to the last the updates (scan.stat.zdry re-added, treatment of GEOCENTR, etc.)
%   09 Feb 2017 by D. Landskron: application of GRAD corrected
%   13 Feb 2017 by A. Hellerschmied: bug-fixes ("special_handling_tag" and calculation of station long/lat for obs. epochs)
%   20 Feb 2017 by H. Krasna: Antenna axis offset altitude correction added
%   23 Feb 2017 by A. Hellerschmied: flagmess is no global variable any more
%   10 Mar 2017 by D. Landskron: command window output for GPT3 backup corrected
%   14 Mar 2017 by A. Hellerschmied: Session name now taken from "parameter.session_name"
%   06 Jun 2017 by A. Hellerschmied: Corrected calculation of partial derivatives for PWL satellite position offsets
%   22 Jun 2017 by A. Hellerschmied: "psou" was not stored in scan structure
%   13 Sep 2017 by D. Landskron: 'tropSource' shifted into 'vie_init' 
%   11 May 2018 by D. Landskron: bug with usage of raytr-files corrected
%   05 Jul 2018 by D. Landskron: vm1 renamed to vmf1 and VMF3 added to the troposphere models 
%   27 Jul 2019 by D. Landskron: zwet parameter added to scan structure
%   15 Jan 2020 by M. Mikschi: gravitational deformation correction added
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%  NOTATION:
%               UPPER CASE ........ variable vector for all scans
%               lower case ........ one value (vector) for one scan
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [antenna, sources, scan, parameter] = vie_mod(antenna, sources, scan, parameter)


% ##### Options #####

% Set flag = true to save results in variable "erg" ("erg" = "Ergebnisse" = "results") for debugging, etc.
flag_save_results = false;
flag_fix_sat_pos_to_stat1 = 0;

% Check special handling tag in parameter structure:
% => initialize it, if it does not exist
if isfield(parameter.vie_mod, 'special_handling_tag')
    if strcmp(parameter.vie_mod.special_handling_tag, 'geocentricAprioriDelay')
        flag_save_results = true;
    end
else
    parameter.vie_mod.special_handling_tag = 'none';
end
    
if flag_save_results
	i_result = 0;
end

% Select dealy model for quasars:
%   1 = Consensus model (standard in VLBI)
%   2 = Sovers / influence of planets not included yet
%   3 = SSB acceleration (calculation of the time delay following the consensus model + Titov 2011)
del_mod_q = 1;

% Select delay model for space crafts:
%   1 = Light time equation (geocentric model)
del_mod_s = 1;

% Options for iteration of sat. velocity and position:
ddt_threshold            = 1e-12;   % Iteration ends, if the change in "dt" is less than this value [sec]
max_iterations          = 10;       % Max. number of iterations


disp('---------------------------------------------------------------')
disp('|                     Welcome to VIE_MOD                      |')
disp('---------------------------------------------------------------')
fprintf('\n')

opt=parameter.lsmopt;

% ##### Check input data #####
% If the observation type file is missing all scans are supposed to be quasar scans!
if ~isfield(scan, 'obs_type')
    [scan.obs_type] = deal('q');
end

% ��������������������������
% �  0 Prepare data                                                      �
% ��������������������������

constants
global c omega
global gms gme gmm

% Get session name:
session = parameter.session_name;

if parameter.vie_mod.write_jet==1
    fidjet = fopen(parameter.vie_init.jetfilnam,'w');
    fidjetuv = fopen(parameter.vie_init.jetfilnamuv,'w');
    fidjetjb = fopen(parameter.vie_init.jetfilnamjb,'w');
    fprintf('writing jetangles to file\n');
end

% ##### Evaluate the GPT3 to get meteorological data #####
% => Update the scan and antenna strucure with met. data (and mapping function parameters from GPT3)
[antenna, scan] = calc_met_data(parameter, antenna, scan);


% ##### CW Output in case met. data from backup models was used #####
% Temp.:
if strcmp(parameter.vie_init.tp,'in situ')  % this if statement is only used in terms of speed, calculation would be the same without
    out_str = {};
    for i_stat = 1 : length(antenna)
        if antenna(i_stat).gpt3temp
            out_str{length(out_str) + 1} = antenna(i_stat).name;
        end
    end
    if ~isempty(out_str)
        fprintf(1, ' => Temperature values from GPT3 (backup) were used for:\n');
        for i_tmpt = 1 : length(out_str)
            fprintf(1, '     - %s\n', out_str{i_tmpt});
        end
    end
end
% Pres.:
if strcmp(parameter.vie_init.zhd,'in situ')
    out_str = {};
    for i_stat = 1 : length(antenna)
        if antenna(i_stat).gpt3pres
            out_str{length(out_str) + 1} = antenna(i_stat).name;
        end
    end
    if ~isempty(out_str)
        fprintf(1, ' => Pressure values from GPT3 (backup) were used for:\n');
        for i_tmpt = 1 : length(out_str)
            fprintf(1, '     - %s\n', out_str{i_tmpt});
        end
    end
end
% e.:
if strcmp(parameter.vie_init.zwd,'in situ')
    out_str = {};
    for i_stat = 1 : length(antenna)
        if antenna(i_stat).gpt3e
            out_str{length(out_str) + 1} = antenna(i_stat).name;
        end
    end
    if ~isempty(out_str)
        fprintf(1, ' => Rel. humidity values (e) from GPT3 (backup) were used for:\n');
        for i_tmpt = 1 : length(out_str)
            fprintf(1, '     - %s\n', out_str{i_tmpt});
        end
    end
end

% Get total number of scans
number_of_all_scans  = length(scan);

% Get RaDec coordinates and names of observed quasars:
RA2000      = zeros(number_of_all_scans,1);
DE2000      = zeros(number_of_all_scans,1);
sourceNames = cell(number_of_all_scans,1);

obs_type_q_ind = strcmp({scan.obs_type}, 'q');
obs_type_s_ind = strcmp({scan.obs_type}, 's');

if sum(obs_type_q_ind) ~= 0
    RA2000(obs_type_q_ind)          = deal([sources.q([scan(obs_type_q_ind).iso]).ra2000]);
    DE2000(obs_type_q_ind)          = deal([sources.q([scan(obs_type_q_ind).iso]).de2000]);

    if strncmpi('icrf3',parameter.vie_init.crf(2),5) && (del_mod_q ~= 3)
        [DE2000, RA2000] = correct_GA(DE2000,RA2000,mean([scan(:).mjd]));
        fprintf(1, 'ICRF3 is used --> GA will be corrected to 2015 using 5.8 muas/year\n');
    end
    sourceNames(obs_type_q_ind)     = deal({sources.q([scan(obs_type_q_ind).iso]).name});
end
if sum(obs_type_s_ind) ~= 0
    sourceNames(obs_type_s_ind)     = deal({sources.s([scan(obs_type_s_ind).iso]).name});
end

% % Number of scans (of q and s type):
% number_of_q_scans   = sum(obs_type_q_ind);
% number_of_s_scans   = sum(obs_type_s_ind);

% Prepare time vectors (observation epochs):
MJD         = ([scan.mjd])';                % time of observations MJD UT
LEAP        = tai_utc(MJD);                 % get difference TAI-UTC [s]
TT          = MJD + (32.184 + LEAP)./86400; % [MJD TT]



% �������������������������� 
% �  1. EARTH ORIENTATION  �
% ��������������������������

% -------------------------------------------------
%   READING EARTH ORIENTATION PARAMETERS
% -------------------------------------------------

% Get EOP values for observation epochs:
[MJDeop, XPeop, YPeop, UT1eop, dXeop, dYeop] = load_eop(MJD, parameter); % [d,rad,sec]

% store a priori values at EOP timesteps
parameter.eop.mjd =  MJDeop;
parameter.eop.xp  =  rad2as(XPeop); % [as]
parameter.eop.yp  =  rad2as(YPeop); % [as]
parameter.eop.ut1 =         UT1eop; % [sec]
parameter.eop.dX  =  rad2as(dXeop); % [as]
parameter.eop.dY  =  rad2as(dYeop); % [as]

% Get EOP values for orbit data epochs:
if ~isempty(sources.s)
    % Get time epochs of orbit data in sources structure:
    MJD_s     = [sources.s(1).mjd];                     % [MJD UTC]
    LEAP_s    = tai_utc(MJD_s);                         % get difference TAI-UTC [s]
    TT_s      = MJD_s + (32.184 + LEAP_s)./86400;       % [MJD TT]
    [MJDeop_s, XPeop_s, YPeop_s, UT1eop_s, dXeop_s, dYeop_s] = load_eop(MJD_s, parameter); % [d,rad,sec]
    number_of_orbit_epochs = length(MJD_s);
end


% -------------------------------
%  INTERPOLATE EOP
% -------------------------------

% ##### Interpolate EOP values for observation epochs: #####
disp('Interpolate EOP values for observation epochs:')

% subtraction of tidal variations (Defraigne and Smits) in dUT1 before interpolation:
if parameter.vie_mod.tidalUT == 1
%     disp('remove tidal UT')
    taiut    = tai_utc(MJDeop);
    MJDTTeop = MJDeop + (32.184 + taiut)/86400;
    %         UT1corr  = tver2000(MJDTTeop);  % [sec]
    if parameter.vie_mod.tidalUT35 ==1
        par35=1;
    else
        par35=2;
    end
    UT1corr  = rg_zont2(MJDTTeop, par35);  % [sec]
    UT1eop   = UT1eop - UT1corr;
end

% Interpolation:
% Linear
if parameter.vie_mod.linear == 1
    disp('linear interpolation of EOP')
    parameter.eop.interp = 'linear';
    if parameter.vie_mod.linear48h
        disp('special linear interpolation for IVS SINEX submission with 48h EOP interval')
        % find index of MJDeop that lies within the 48 hour estimation interval
        midof48h = find(MJDeop==floor(min(MJD))+1);
        % remove this midnight point from a priori time series
        MJDeop(midof48h) = []; UT1eop(midof48h) = [];
        XPeop(midof48h) = []; YPeop(midof48h) = []; 
        dXeop(midof48h) = []; dYeop(midof48h) = [];

        % store a priori values at EOP timesteps without the 24h midnight
        parameter.eop.mjd(midof48h) =  [];
        parameter.eop.xp(midof48h) =  [];
        parameter.eop.yp(midof48h) =  [];
        parameter.eop.ut1(midof48h) =  [];
        parameter.eop.dX(midof48h) =  [];
        parameter.eop.dY(midof48h) =  [];
    end
    % A priori EOP values are determined with linear interpolation between
    % the value of midnight before and after observation time.
    % for a session from 18:00 to 18:00 this means, that there are 2 a
    % priori lines and a break at midnight
    % no a priori values are stored in the parameter file!!!
    DUT1  = interp1(MJDeop,UT1eop,MJD,'linear','extrap');
    XP    = interp1(MJDeop, XPeop,MJD,'linear','extrap');
    YP    = interp1(MJDeop, YPeop,MJD,'linear','extrap');
    DX    = interp1(MJDeop, dXeop,MJD,'linear','extrap');
    DY    = interp1(MJDeop, dYeop,MJD,'linear','extrap');

% Lagragne
else % linear = 0
    disp('Lagrange interpolation of EOP')
    parameter.eop.interp = 'lagrange';
    % subtraction of tidal variations (Defraigne and Smits) in dUT1 before interpolation
    % interpolate EOP for time of observation
    DUT1    = lagint4v(MJDeop,UT1eop,MJD);
    XP      = lagint4v(MJDeop, XPeop,MJD);
    YP      = lagint4v(MJDeop, YPeop,MJD);
    DX      = lagint4v(MJDeop, dXeop,MJD);
    DY      = lagint4v(MJDeop, dYeop,MJD);
end

% re-add tidal variation in dUT1 after interpolation:
if parameter.vie_mod.tidalUT == 1
%     disp('re-add tidal UT')
    corrUT1 = rg_zont2(TT,par35);  % [sec]
    DUT1    = DUT1 + corrUT1;   % [sec]
end

% Add high frequency EOP
[CORX,CORY,CORUT,parameter] = eophf(TT, parameter);
DUT1 = DUT1 + CORUT;
XP   = XP + CORX;
YP   = YP + CORY;

disp('...done.')


% ##### Interpolate EOP values for orbit data epochs: #####
if ~isempty(sources.s)
    disp('Interpolate EOP values for satellite orbit epochs:')
    
    % subtraction of tidal variations (Defraigne and Smits) in dUT1 before interpolation:
    if parameter.vie_mod.tidalUT == 1
        disp('remove tidal UT')
        taiut    = tai_utc(MJDeop_s);
        MJDTTeop_s = MJDeop_s + (32.184 + taiut)/86400;
        %         UT1corr  = tver2000(MJDTTeop);  % [sec]
        if parameter.vie_mod.tidalUT35 ==1
            par35 = 1;
        else
            par35 = 2;
        end
        UT1corr_s  = rg_zont2(MJDTTeop_s, par35);  % [sec]
        UT1eop_s   = UT1eop_s - UT1corr_s;
    end
    
    % Interpolation:
    % Linear
    if parameter.vie_mod.linear == 1
        disp('linear interpolation of EOP')
        parameter.eop.interp = 'linear';
        % A priori EOP values are determined with linear interpolation between
        % the value of midnight before and after observation time.
        % for a session from 18:00 to 18:00 this means, that there are 2 a
        % priori lines and a break at midnight
        % no a priori values are stored in the parameter file!!!
        DUT1_s  = interp1(MJDeop_s,UT1eop_s, MJD_s,'linear', 'extrap');
        XP_s    = interp1(MJDeop_s, XPeop_s, MJD_s,'linear', 'extrap');
        YP_s    = interp1(MJDeop_s, YPeop_s, MJD_s,'linear', 'extrap');
        DX_s    = interp1(MJDeop_s, dXeop_s, MJD_s,'linear', 'extrap');
        DY_s    = interp1(MJDeop_s, dYeop_s, MJD_s,'linear', 'extrap');
        % Lagragne
    else % linear = 0
        disp('Lagrange interpolation of EOP')
        parameter.eop.interp = 'lagrange';
        % subtraction of tidal variations (Defraigne and Smits) in dUT1 before interpolation
        % interpolate EOP for time of observation
        DUT1_s    = lagint4v(MJDeop_s,UT1eop_s, MJD_s);
        XP_s      = lagint4v(MJDeop_s, XPeop_s, MJD_s);
        YP_s      = lagint4v(MJDeop_s, YPeop_s, MJD_s);
        DX_s      = lagint4v(MJDeop_s, dXeop_s, MJD_s);
        DY_s      = lagint4v(MJDeop_s, dYeop_s, MJD_s);
    end
    
    % re-add tidal variation in dUT1 after interpolation:
    if parameter.vie_mod.tidalUT == 1
        disp('re-add tidal UT')
        corrUT1_S = rg_zont2(TT_s, par35);  % [sec]
        DUT1_s    = DUT1_s + corrUT1_S;   % [sec]
    end
    
    % Add high frequency EOP
    [CORX_s, CORY_s, CORUT_s, parameter] = eophf(TT_s, parameter);
    DUT1_s = DUT1_s + CORUT_s;
    XP_s   = XP_s + CORX_s;
    YP_s   = YP_s + CORY_s;
    
    disp('...done.')
    
end % if ~isempty(sources.s)


% ##### if ray-tracing files are to be used #####
raytr_used = strcmp(parameter.vie_init.tropSource.name,'raytr');

% if ray-traced delays are used, then set the remaining parameters to VMF3 which can then be used as default, in case a certain line is missing in the .trp-file
if raytr_used
    
    parameter.vie_init.zhd = 'in situ';
    parameter.vie_init.zwd = 'no';
    parameter.vie_mod.mfh = 'vmf3';
    parameter.vie_mod.mfw = 'vmf3';
    parameter.vie_mod.apgm_h = 'no';
    parameter.vie_mod.apgm_w = 'no';
    
end



% -----------------------------------------------------------------------
% barycentric coordinates & velocities of celestial bodies (EPHEMERIDES)
% -----------------------------------------------------------------------
ephnum = parameter.vie_mod.eph;
ephnam=([ephnum,'_',session,'.mat']);
ephok=0;
if exist(['../EPHEM/',ephnam],'file')
    load(['../EPHEM/',ephnam],'ephem')
    disp('load existing ephemerides ...')
    if size(ephem.time) == size(TT)
        if ephem.time == TT
            ephok=1;
        end
    end
end
if ephok==0
    disp(strcat('reading_ ' , ephnum , ' ephemerides'))
    ephem = load_eph(TT,ephnum);
    save(['../EPHEM/',ephnam],'ephem')
end

% ------------------------------------------------------------------------
% Transformation matrix Q TRS --> CRS following Conv. Chapt. 5
%    + partial derivatives dQ/dEOP
% ------------------------------------------------------------------------
nutmod = parameter.vie_mod.nutmod;

% special call of trs2crs only needed for tidal ERP variations partials
if opt.est_erpsp==1
    [T2C,DQDXP,DQDYP,DQDUT,DQDDX,DQDDY,XNUT,YNUT,ERA,PNN,RN,WN,R3SS,R2XP,DR2XP,R1YP,DR1YP,DRDUTN,COSAG1,SINAG1,COSAG2,SINAG2] = trs2crs_tidERP(MJD,XP,YP,DUT1,DX,DY,nutmod);
    % special call of trs2crs only needed for FCN partials
elseif opt.est_FCNnut == 1 || opt.est_FCNnutAmpl == 1
    [T2C,DQDXP,DQDYP,DQDUT,DQDDX,DQDDY,XNUT,YNUT,ERA,DQDFCN,DQDFCNAC,DQDFCNAS] = trs2crs_fcn(MJD,XP,YP,DUT1,DX,DY,nutmod);
else
    [T2C,DQDXP,DQDYP,DQDUT,DQDDX,DQDDY,XNUT,YNUT,ERA] = trs2crs(MJD,XP,YP,DUT1,DX,DY,nutmod);
end

% Get trafo. matrices (T2C_s) to convert orbit data from TRS (e.g. from SP3 files) to CRS:
if ~isempty(sources.s)
    [T2C_s, DQDXP_s, DQDYP_s, DQDUT_s, DQDDX_s, DQDDY_s, XNUT_s, YNUT_s, ERA_s] = trs2crs(MJD_s, XP_s, YP_s, DUT1_s, DX_s, DY_s, nutmod);
end % if ~isempty(sources.s)


% --------------------------------------------------------------
%  source vector and derivatives in the catalogue system
% --------------------------------------------------------------

% #### Quasars: ####
[RQ,DRQDRA,DRQDDE] = sourcevec(DE2000,RA2000); % Output: [source vector, partial derivative w.r.t. ra, partial derivative w.r.t. de]

% #### Space Crafts: ####
% Source vectors:
if ~isempty(sources.s)
    % loop over space crafts:
    for i_sc = 1 : length(sources.s)
        xyz_crf_tmp     = zeros(number_of_orbit_epochs, 3);
        if  sources.s(i_sc).flag_v_trf
            v_xyz_crf_tmp   = zeros(number_of_orbit_epochs, 3);
        end
        % Loop over all oribit pos. epochs:
        for i_orbit_epoch = 1 : number_of_orbit_epochs
            % Position:
            xyz_crf_tmp(i_orbit_epoch, :) = (T2C_s(:, :, i_orbit_epoch) * [sources.s(i_sc).x_trf(i_orbit_epoch); sources.s(i_sc).y_trf(i_orbit_epoch); sources.s(i_sc).z_trf(i_orbit_epoch)])';
            % Velocity:
            if  sources.s(i_sc).flag_v_trf
                v_xyz_crf_tmp(i_orbit_epoch, :) = (T2C_s(:, :, i_orbit_epoch) * ([sources.s(i_sc).vx_trf(i_orbit_epoch); sources.s(i_sc).vy_trf(i_orbit_epoch); sources.s(i_sc).vz_trf(i_orbit_epoch)] + cross([0; 0; omega], [sources.s(i_sc).x_trf(i_orbit_epoch); sources.s(i_sc).y_trf(i_orbit_epoch); sources.s(i_sc).z_trf(i_orbit_epoch)])))';
            end
        end
        sources.s(i_sc).x_crf = xyz_crf_tmp(:, 1);
        sources.s(i_sc).y_crf = xyz_crf_tmp(:, 2);
        sources.s(i_sc).z_crf = xyz_crf_tmp(:, 3);
        if  sources.s(i_sc).flag_v_trf
            sources.s(i_sc).vx_crf = v_xyz_crf_tmp(:, 1);
            sources.s(i_sc).vy_crf = v_xyz_crf_tmp(:, 2);
            sources.s(i_sc).vz_crf = v_xyz_crf_tmp(:, 3);
        end
    end
end % if ~isempty(sources.s)
% Derivatives:
% => Calculated later on!


% #### Source structure ####
% choose simulated structure catalogue depending on options
if parameter.vie_mod.ssou==1 || parameter.vie_mod.write_jet==1 % have structure
    % read in catalogue
    filename_sou_cat = strcat('../CRF/SOURCE_STRUCTURE_CAT/',parameter.vie_mod.sou_cat);
    [cat_comp.name,cat_comp.flux,cat_comp.maj,cat_comp.min,cat_comp.angle,cat_comp.dRA,cat_comp.dDec] = textread(filename_sou_cat,'%s%f%f%f%f%f%f','delimiter',',');
end




% ������������������������������������������������������������������������
% �  2. STATION COORDINATES                                              �
% ������������������������������������������������������������������������
zz = zeros(length(antenna),1);
flagmess = struct('cto',zz,'cta',zz,'cnta',zz,'crg',zz,'gia',zz,'axtyp',zz,'thermal',zz,'vmf3',zz,'vmf1',zz,'dao',zz,'ctop',zz,'chl',zz);

% antenna corrections for all stations
time_lim = [min(MJD) max(MJD)];     % scans might not be in straight order
[antenna, flagmess] = corr_ant(time_lim, scan(1).tim, antenna, parameter, flagmess);

% antenna positions and velocities
ANT = [[antenna.x]', [antenna.y]', [antenna.z]'];
VEL = [[antenna.vx]', [antenna.vy]', [antenna.vz]'];


% a priori gradients

%year = mjd2date(MJD(1));   % there are no multi-year sessions
[year, month, day, ~, ~, ~] = mjd2date(MJD(1));

% Preallocate a priori gradients for each antenna (north and east gradient)
% - initialize with zeros
a_ngr_h = zeros(length(antenna), 1);
a_egr_h = zeros(length(antenna), 1);
a_ngr_w = zeros(length(antenna), 1);
a_egr_w = zeros(length(antenna), 1);
a_ngr = zeros(length(antenna), 1);
a_egr = zeros(length(antenna), 1);


% define relevant GRAD data, if needed
if strcmpi(parameter.vie_mod.apgm_h,'grad')   ||   strcmpi(parameter.vie_mod.apgm_w,'grad')
    
    % read respective file
    gradPath = '../TRP/GRAD/';
    gradFile = [gradPath, 'y', num2str(year), '.grad_r'];
    fidGrad = fopen(gradFile, 'r');
    if fidGrad ~= -1
        gradData = textscan(fidGrad, '%8s  %8.2f  %6.3f %6.3f %6.3f %6.3f', 'CommentStyle', '#','delimiter', '||');
        fclose(fidGrad);
        % the following is not necessary
        %         % if the date is close to a turn of the year, then read also the following/preceding year
        %         if month==1 && day < 2
        %             gradFile = [gradPath, 'y', num2str(year-1), '.grad_r'];
        %             fidGrad = fopen(gradFile, 'r');
        %             gradData_add = textscan(fidGrad, '%8s  %8.2f  %6.3f %6.3f %6.3f %6.3f', 'CommentStyle', '#','delimiter', '||');
        %             fclose(fidGrad);
        %             gradData_all = [gradData;gradData_add];
        %         end
        %         if month==12 && day > 30
        %             gradFile = [gradPath, 'y', num2str(year+1), '.grad_r'];
        %             fidGrad = fopen(gradFile, 'r');
        %             gradData_add = textscan(fidGrad, '%8s  %8.2f  %6.3f %6.3f %6.3f %6.3f', 'CommentStyle', '#','delimiter', '||');
        %             fclose(fidGrad);
        %             gradData_all = [gradData;gradData_add];
        %         end
        %         % if another year was appended, data needs to be sorted
        %         if exist('gradData_all','var')
        %             for i = 1:size(gradData_all,2)
        %                 gradData{i} = cat(1, gradData_all{:,i});
        %             end
        %             [~,ind] = sort(gradData{2});
        %             for i = 1:size(gradData,2)
        %                 gradData{i} = gradData{i}(ind);
        %             end
        %
        %         end
        %
        %         % get only wanted lines (i.e. lines with only session-included stations)
        %         wantedLines = ismember(gradData{1} , {antenna.name});
        %         for k=1:size(gradData,2)
        %             gradData{k}=gradData{k}(wantedLines);
        %         end
        %
        %         % find and exclude duplicates (same station and same epoch)
        %         if exist('gradData_all','var')
        %             [~,ind_withoutDupl,~] = unique(strcat(gradData{1},num2str(gradData{2})));
        %             for i = 1:length(gradData)
        %                 gradData{i} = gradData{i}(sort(ind_withoutDupl));
        %             end
        %         end
    else
        fprintf('The required GRAD file is not available, no a priori gradients used instead\n');
        parameter.vie_mod.apgm_h = 'no';
        parameter.vie_mod.apgm_w = 'no';
    end
end
    

% Elliptical coordinates of antennas
[PHI,LAM,H_ELL] = xyz2ell(ANT);


% acceleration of SSB - Galactocentric aberration
% GA value:
GA=5.8; %[uas/y]
RAGC = (17+45/60+40/3600)/12*pi; %[rad] RA of Galactic center: 17h45min40sec
DeGC = (-29-28/3600)/180*pi;   %[rad]  De of Galactic center: -29deg 00min 28sec

% num of sec in one year
GA_T = 60*60*24*365.25;
% conversion rad -> as
GA_k=180/pi*60*60 *GA_T; % as*s
accMag = GA* 1e-6 *c  / GA_k ; % m/s^2 magnitude of acceleration (~ 2.44e-10 m/s^2)

acc=[accMag*cos(RAGC)*cos(DeGC)
     accMag*sin(RAGC)*cos(DeGC)
     accMag*sin(DeGC)]; %m/sec^2




% *************************
%  loop over scans
% *************************
disp('station corrections')

% do that before the loop - otherwise very likely pretty slow!
cpsd_all = cPostSeismDeform(MJD,antenna); % [3 x nScans x nStat] matrix / meters!

if strcmpi(parameter.vie_init.zhd,'gpt3')   ||   strcmpi(parameter.vie_init.zwd,'gpt3')   ||   strcmpi(parameter.vie_mod.mfh,'gpt3')   ||   strcmpi(parameter.vie_mod.mfw,'gpt3')   ||   strcmpi(parameter.vie_mod.apgm_h,'gpt3')   ||  strcmpi(parameter.vie_mod.apgm_w,'gpt3')
    cell_grid_GPT3 = gpt3_5_fast_readGrid;
end

% GRAVITATIONAL DEFORMATION
% check if the antenna struct has a field for gravitational deformation. It
% might have been created with an older version of vie_init. If not warn
% the user.
if parameter.vie_mod.gravDef == 1 && ~isfield(antenna, 'gravdef')
    fprintf(['WARNING: Correction for gravitational antenna deformation is \n'...
          'turned on but the used antenna struct has no gravdef field!\n'])
end
% check if any of the stations in the session has available gravdef data.
% The user might use an old superstation-file without this information.
if parameter.vie_mod.gravDef == 1 && isfield(antenna, 'gravdef') && ...
        isempty([antenna.gravdef])
    fprintf(['WARNING: Correction for gravitational antenna deformation \n'...
          'is turned on but none of the stations in the session has \n'...
          ' gravdef information in the \n antenna struct!\n'])
end

for isc = 1:number_of_all_scans
    
    % running variables for active scan
    mjd   = MJD(isc);
    tim   = scan(isc).tim;
    leap  = LEAP(isc);
    t2c   = T2C(:,:,isc);
    dQdxp = DQDXP(:,:,isc);
    dQdyp = DQDYP(:,:,isc);
    dQdut = DQDUT(:,:,isc);
    dQddX = DQDDX(:,:,isc);
    dQddY = DQDDY(:,:,isc);
    
    sec_of_day = scan(isc).tim(4)*3600 + scan(isc).tim(5)*60 + scan(isc).tim(6);
    
    % tidal ERP variations
    if opt.est_erpsp==1
        PNn     = PNN(:,:,isc);
        Rall    = RN(:,:,isc);
        Wn      = WN(:,:,isc);
        R3ss    = R3SS(:,:,isc);
        R2xp    = R2XP(:,:,isc);
        dR2xp   = DR2XP(:,:,isc);
        R1yp    = R1YP(:,:,isc);
        dR1yp   = DR1YP(:,:,isc);
        dRdutn  = DRDUTN(:,:,isc);
        cosag2  = COSAG2(:,isc)';
        sinag2  = SINAG2(:,isc)';
        cosag   = [COSAG1(:,isc)' COSAG2(:,isc)'];
        sinag   = [SINAG1(:,isc)' SINAG2(:,isc)'];
    end
    
    % FCN partials from nutation
    if opt.est_FCNnut == 1
        dQdfcn = DQDFCN(:,:,isc);
    end
    if opt.est_FCNnutAmpl == 1
        dQdfcnAc = DQDFCNAC(:,:,isc);
        dQdfcnAs =DQDFCNAS(:,:,isc);
    end
    
    % Source:
    % RQ, DRQDRA, DRQDDE only contain valid values for quasar observations.
    % => For spacecraft obs. the values are = 0!
    rq      = RQ(isc,:);         % Source vector
    drqdra  = DRQDRA(isc,:);     % Partial derivative w.r.t. Ra
    drqdde  = DRQDDE(isc,:);     % Partial derivative w.r.t. Dec
    
    % Pole coordinates:
    xp      = XP(isc);
    yp      = YP(isc);
    
    % Ephemerids:
    earth   = ephem.earth(isc).xbar;
    vearth  = ephem.earth(isc).vbar;
    moon    =  ephem.moon(isc).xgeo;
    sun     =   ephem.sun(isc).xgeo;
    
    % Compute frequency and phase of 342 tidal constituents, which are
    % recommended for implementation of ocean tidal loading in IERS
    % Conventions 2010.
    [cto_F, cto_P, cto_TAMP, cto_IDD1] = libiers_tdfrph_call(mjd, leap);

    % reference epoch for SSB acceleration
    refep_accSSB = 57023; % 2015.0
    delt_accSSB  = (mjd - refep_accSSB) * 86400; % time since ref epoch in sec
 
    % ***********************
    %  loop over stations in current scan
    % ***********************
    for i_stat = 1 : length(scan(isc).stat)
        
        % Choose only active stations in the scan:
        if ~isempty(scan(isc).stat(i_stat).temp)
            ant = ANT(i_stat,:);
            vel = VEL(i_stat,:);
            
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
                if isempty(antenna(i_stat).cto) == 0
                    cto_data = antenna(i_stat).cto;
                    %[cto] = ctocean(mjd,leap,ant(ist,:),cto_data);
                    [cto_rsw] = libiers_hardisp(mjd, leap, cto_data, cto_F,cto_P, cto_TAMP, cto_IDD1); % [Radial, South, West];  IERS Conv. 2010
                    cto = ren2xyz([cto_rsw(1) -cto_rsw(3) -cto_rsw(2)] , PHI(i_stat), LAM(i_stat));
                end
            end
            
            % Tidal Atmosphere Loading
            cta = [0 0 0];
            if parameter.vie_mod.cta == 1
                if isempty(antenna(i_stat).cta) == 0
                    S12_data = antenna(i_stat).cta;
                    S12_data = S12_data(1:12);         % for "ATIDE_LEONID.mat"
                    [cta] = ctatmos12(mjd, ant, S12_data);
                end
            end
            
            % Non-tidal Atmosphere Loading
            cnta = [0 0 0];
            if parameter.vie_mod.cnta == 1
                if isempty(antenna(i_stat).cnta_dx) == 0
                    cnta_data = antenna(i_stat).cnta_dx;
                    
                    cnta(:,1) = call_spline_4(cnta_data(:,1), cnta_data(:,2), mjd);
                    cnta(:,2) = call_spline_4(cnta_data(:,1), cnta_data(:,3), mjd);
                    cnta(:,3) = call_spline_4(cnta_data(:,1), cnta_data(:,4), mjd);
                end
            end
            
            crg  = [0 0 0]; %xyz [m]
            cgia = [0 0 0]; %xyz [m]
            
            % Atmosphere loading from regression coefficients
            pres = scan(isc).stat(i_stat).pres;  % total pressure [hPa]
            pres0 = antenna(i_stat).crg(1); %[hPa]
            if parameter.vie_mod.crg == 1
                rg = antenna(i_stat).crg(2); %[m/hPa]
                crg_h = (pres - pres0) * rg;
                crg = ren2xyz([crg_h 0 0], PHI(i_stat), LAM(i_stat));
            end
            pcrg_h      = (pres-pres0); %[hPa]
            pcrg_xyz    = ren2xyz([pcrg_h 0 0],PHI(i_stat),LAM(i_stat));  %[hPa]
            
            
            % GIA uplift rates
            if parameter.vie_mod.gia == 1
                vgia = antenna(i_stat).gia_dvx(2:4); % m/y
                mjd0_gia = antenna(i_stat).gia_dvx(1); %mjd
                
                cgia = vgia /365.25 * (mjd - mjd0_gia); %m
            end
            
            
            
            % Hydrology loading
            chl = [0 0 0];
            if parameter.vie_mod.chl == 1
                if isempty(antenna(i_stat).chl_dx) == 0
                    chl_data = antenna(i_stat).chl_dx;
                    chl(:,1) = call_spline_4(chl_data(:,1), chl_data(:,2), mjd);
                    chl(:,2) = call_spline_4(chl_data(:,1), chl_data(:,3), mjd);
                    chl(:,3) = call_spline_4(chl_data(:,1), chl_data(:,4), mjd);
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
                if isempty(antenna(i_stat).opl) == 0
                    opl = antenna(i_stat).opl;
                    ctpm = parameter.vie_mod.ctpm; % mean: linear or cubic
                    [ctop, flgm_ctp] = ctoceanpole(tim, ant, xp, yp, opl, ctpm);
                end
            end
            
            % Seasonal variations of station positions
            % (zero a priori; only partial derivatives w.r.t. amplitudes are build up)
            pAcr_xyz = [0 0 0]; pAce_xyz = [0 0 0]; pAcn_xyz = [0 0 0];
            pAsr_xyz = [0 0 0]; pAse_xyz = [0 0 0]; pAsn_xyz = [0 0 0];
            if opt.est_stsespos
                [pAcr_xyz, pAce_xyz, pAcn_xyz, pAsr_xyz, pAse_xyz, pAsn_xyz]= stseasonal(mjd,PHI(i_stat),LAM(i_stat));
            end
            
            % post-seismic deformation (ITRF2014 at least)
            cpsd = cpsd_all(:, isc,i_stat)'; % [3 x nScans x nStat] matrix / meters!
            
            % Displacement of reference markers on the crust
            % => For the GEOCENTR station displacements are set to zero!
            if strcmp(antenna(i_stat).name,'GEOCENTR')
                cposit = 0;
            else
                cposit = cts + cto + cta + cnta + crg + cgia + ctp + ctop + chl + cpsd;
            end
            
            % time since reference epoch [years]
            tep = (mjd - antenna(i_stat).epoch)/365.25;
            
            % Real antenna position
            antv     = ant + tep * vel;
            ant_trs  = antv + cposit + antenna(i_stat).c_ecc;
            ant_crs  = t2c * ant_trs'; % transformation to the cel. system
            
            % store results per scan & station
            scan(isc).stat(i_stat).x      = ant_trs; % corrected station position in TRS [x,y,z]
            scan(isc).stat(i_stat).xcrs   = ant_crs; % corrected station position in CRS [x,y,z]
            
            scan(isc).stat(i_stat).pAcr_xyz = pAcr_xyz; % [-] partials for cosine Ampl. - seasonal stat. position variations
            scan(isc).stat(i_stat).pAce_xyz = pAce_xyz; % [-]
            scan(isc).stat(i_stat).pAcn_xyz = pAcn_xyz; % [-]
            scan(isc).stat(i_stat).pAsr_xyz = pAsr_xyz; % [-] partials for sine Ampl. - seasonal stat. position variations
            scan(isc).stat(i_stat).pAse_xyz = pAse_xyz; % [-]
            scan(isc).stat(i_stat).pAsn_xyz = pAsn_xyz; % [-]
            
            scan(isc).stat(i_stat).phpole = phpole;     % [cm] partials for pole tide Love numbers
            scan(isc).stat(i_stat).plpole = plpole;     % [cm]
            
            scan(isc).stat(i_stat).prg = pcrg_xyz;      % [hPa] % partials for APL RG
            
            scan(isc).stat(i_stat).pLove    = pdxyz_h;  % [cm] partials for Love numbers
            scan(isc).stat(i_stat).pShida   = pdxyz_l;  % [cm] partials for Shida numbers
            scan(isc).stat(i_stat).pFCN     = pdxyz_FCN;% [cm] partials for FCN
            
            
            % -------------------------
            %  tropospheric parameters (using ray-traced delays)
            % -------------------------
            
            if raytr_used   % if ray-tracing files are used
                if ~exist('firstRaytrRun', 'var')   % read the .trp-file only in the first run
                    disp('use ray-tracing files')
                    firstRaytrRun = 1;
                    [raytrdata, raytrFileFoundLog] = load_trpfile(parameter, session);
                end
                scan = get_trpdel(raytrdata, scan, isc, i_stat, antenna, raytrFileFoundLog, sourceNames);
            end
            
            aht = 0; awt = 0; zhdt = 0; zwdt = 0;
            if strcmpi(parameter.vie_init.tropSource.name,'indModeling')   ||   isempty(scan(isc).stat(i_stat).trop)
                
                if strcmpi(parameter.vie_init.zhd,'vmf3')
                    tmjd = antenna(i_stat).vmf3(:,1);
                    zhd  = antenna(i_stat).vmf3(:,4);

                    zhdt = call_spline_4(tmjd,zhd,mjd);
                end
                if strcmpi(parameter.vie_init.zhd,'vmf1')
                    tmjd = antenna(i_stat).vmf1(:,1);
                    zhd  = antenna(i_stat).vmf1(:,4);

                    zhdt = call_spline_4(tmjd,zhd,mjd);
                end
                
                if strcmpi(parameter.vie_init.zwd,'vmf3')
                    tmjd = antenna(i_stat).vmf3(:,1);
                    zwd  = antenna(i_stat).vmf3(:,5);

                    zwdt = call_spline_4(tmjd,zwd,mjd);
                end
                if strcmpi(parameter.vie_init.zwd,'vmf1')
                    tmjd = antenna(i_stat).vmf1(:,1);
                    zwd  = antenna(i_stat).vmf1(:,5);

                    zwdt = call_spline_4(tmjd,zwd,mjd);
                end

                if strcmpi(parameter.vie_mod.mfh,'vmf3')
                    if isempty(antenna(i_stat).vmf3) == 0
                        tmjd = antenna(i_stat).vmf3(:,1);
                        ah   = antenna(i_stat).vmf3(:,2);
                        aw   = antenna(i_stat).vmf3(:,3);
                        
                        aht = call_spline_4(tmjd,ah,mjd);
                        awt = call_spline_4(tmjd,aw,mjd);
                    end
                end
                if strcmpi(parameter.vie_mod.mfh,'vmf1')
                    if isempty(antenna(i_stat).vmf1) == 0
                        tmjd = antenna(i_stat).vmf1(:,1);
                        ah   = antenna(i_stat).vmf1(:,2);
                        aw   = antenna(i_stat).vmf1(:,3);
                        
                        aht = call_spline_4(tmjd,ah,mjd);
                        awt = call_spline_4(tmjd,aw,mjd);
                    end
                end
                
                if strcmpi(parameter.vie_mod.mfw,'vmf3')
                    if isempty(antenna(i_stat).vmf3) == 0
                        tmjd = antenna(i_stat).vmf3(:,1);
                        ah   = antenna(i_stat).vmf3(:,2);
                        aw   = antenna(i_stat).vmf3(:,3);
                        
                        aht = call_spline_4(tmjd,ah,mjd);
                        awt = call_spline_4(tmjd,aw,mjd);
                    end
                end
                if strcmpi(parameter.vie_mod.mfw,'vmf1')
                    if isempty(antenna(i_stat).vmf1) == 0
                        tmjd = antenna(i_stat).vmf1(:,1);
                        ah   = antenna(i_stat).vmf1(:,2);
                        aw   = antenna(i_stat).vmf1(:,3);
                        
                        aht = call_spline_4(tmjd,ah,mjd);
                        awt = call_spline_4(tmjd,aw,mjd);
                    end
                end
                
                
            end % use VMF/GPT3/GMF intead of trp files
            
            % store results per scan & station
            scan(isc).stat(i_stat).aht  = aht;    % interpolated ah
            scan(isc).stat(i_stat).awt  = awt;    % interpolated aw
            scan(isc).stat(i_stat).zhdt = zhdt;   % interpolated zhd
            scan(isc).stat(i_stat).zwdt = zwdt;   % interpolated zwd
            
        end % if ~isempty(scan(isc).stat(i_stat).x)
        
    end % length(scan(isc).stat)
    
    
    % ������������������������������������������������������������������������
    % �  3. COMPUTED DELAY                                                   �
    % ������������������������������������������������������������������������
    
    
    % --------------------------------------
    %  loop over baselines (observations)
    % --------------------------------------
    for iobs = 1:length(scan(isc).obs)
        
        % Loop init.:
        % Partial derivatives of the delay w.r.t. the positon of space crafts (different ref. frames)
        pd_satpos_gcrf  = [];
        pd_satpos_trf   = [];
        pd_satpos_rsw   = [];
        
        % get station IDs of the current baseline:
        stat_1_id  = scan(isc).obs(iobs).i1;
        stat_2_id  = scan(isc).obs(iobs).i2;
        
        % get station information
        stat_1_gcrs = scan(isc).stat(stat_1_id).xcrs;    % station1 position GCRS
        stat_2_gcrs = scan(isc).stat(stat_2_id).xcrs;    % station2 position GCRS
        stat_1_trs  = scan(isc).stat(stat_1_id).x;       % station1 position TRS
        stat_2_trs  = scan(isc).stat(stat_2_id).x;       % station2 position TRS
        
        rqu = rq / norm(rq);                    % unit source vector barycentrum-source, only valid for quasar obs.! For spacecrafts = [1,0,0]
        b_gcrs = stat_2_gcrs - stat_1_gcrs;     % baseline vector GCRS
        b_trs  = stat_2_trs-stat_1_trs;         % baseline vector TRS
        
        % station velocity due to earth rotation
        v1 = [-omega*stat_1_trs(2);omega*stat_1_trs(1);0];  % [TRS]
        v1 = t2c*v1;                                        % [CRS]
        v2 = [-omega*stat_2_trs(2);omega*stat_2_trs(1);0];  % [TRS]
        v2 = t2c*v2;                                        % [CRS]
        
        % delay model   (1) consensus model delmod = 1
        %               (2) Sovers
        %               (3) SSB acceleration 
        %               the accelaration can be also set to zero and
        %               as whole estimated in the global adjustment
        
        
        % ##### Distintguish between source types (quasar/spacecraft) #####
        %   => Different delay models are required!
        switch(scan(isc).obs_type)
            
            % ##### Quasars #####
            case 'q'
                
                % Switch between different delay models:
                switch del_mod_q
                    case 1 % consensus
                        % calculation of the time delay following the consensus model
                        % ------------------------------------------------------------
                        
                        %(1) barycentric station vector (eq. 6)
                        xb1 = earth + stat_1_gcrs;
                        xb2 = earth + stat_2_gcrs;
                        
                        %(2)&(3) differential gravitational delay due to celestial bodies
                        [Tgrav, pGammaSun] = grav_delay(xb1,xb2,vearth,b_gcrs,rqu,ephem,isc,opt);
                        
                        %(4) differential gravitational delay due to the earth
                        Tgrave = 2*gme/c^3*...
                            log((norm(stat_1_gcrs)+rqu*stat_1_gcrs)/(norm(stat_2_gcrs)+rqu*stat_2_gcrs)); % (eq. 1)
                        
                        %(5) total differential gravitational delay (eq. 2)
                        Tgrav = Tgrave + Tgrav;
                        
                        if strcmp(antenna(stat_1_id).name,'GEOCENTR') || strcmp(antenna(stat_2_id).name,'GEOCENTR')
                            Tgrav=0;
                        end
                        
                        %(6) vacuum delay
                        U     = gms/norm(sun); % gravitational potential at the geocenter
                        gamma = 1;
                        fac1  = (rqu*b_gcrs)/c;
                        term1 = 1-(1+gamma)*U/c^2-(norm(vearth))^2/(2*c^2)-(vearth'*v2)/c^2;
                        fac2  = (vearth'*b_gcrs)/c^2;
                        term2 = 1+(rqu*vearth)/(2*c);
                        quot  = 1+(rqu*(vearth+v2))/c;
                        tau   = (Tgrav - fac1*term1 - fac2*term2)/quot;          %(eq. 9)
                        
                        pGammaSun = pGammaSun/quot;
                        
                        %(7) aberrated source vector (eq. 15)
                        k1a = rqu + (vearth+v1)'/c - rqu*((rqu*(vearth+v1))')/c;
                        k2a = rqu + (vearth+v2)'/c - rqu*((rqu*(vearth+v2))')/c;
                        
                    case 2 % sovers / influence of planets not included yet
                        
                        % input
                        bet2s = v2/c;        % velocity st2 GCRS
                        bet   = vearth/c;    % SBB
                        bs    = b_gcrs;           % baseline GCRS
                        k     = -rqu;        % source vector SBB (1,3)
                        
                        %  transform to SBB frame
                        gam  = 1/sqrt(1-bet'*bet);
                        bet2 = (bet2s + (gam-1)*(bet2s'*bet)*bet/norm(bet)^2+gam*bet)/...
                            (gam*(1 + bet2s'*bet)); % v2 SBB
                        b_gcrs    = bs +(gam-1)*(bs'*bet)*bet/norm(bet)^2 -...
                            gam*bet2*(bs'*bet); % baseline SBB
                        
                        % calculate in SBB frame
                        tausbb = (k*b_gcrs)/(1-k*bet2);
                        bt1t2  = b_gcrs + bet2*tausbb;
                        % transform back to GCRS frame
                        tau = gam*tausbb - gam*(bt1t2'*bet);   % [m]
                        tau = tau/c;             % --> [sec]
                        % gravitational correction (Sun)
                        Rr   = ephem.sun(isc).xbar; % position SBB at reception time
                        Vr   = ephem.sun(isc).vbar;
                        trta = norm(earth)/c;
                        Rn   = Rr - trta*Vr/c;
                        % iterative solution for position at time of closest approach
                        Rn = Rr - norm(Rn)*Vr/c;
                        Rn = Rr - norm(Rn)*Vr/c;
                        Rn = Rr - norm(Rn)*Vr/c;
                        Ra = Rn;
                        
                        xb1 = earth + stat_1_gcrs;
                        xb2 = earth + stat_2_gcrs;
                        r1 = xb1 - Ra;
                        r2 = xb2 + bet2*tausbb - Ra;
                        
                        dGp  = 2*gms/c^3*log((norm(r1)+rqu*r1)/(norm(r2)+rqu*r2));
                        U    = gms/(norm(ephem.sun(isc).xgeo)*c^2);
                        dGps = dGp - 2*U*tau;
                        
                        % differential gravitational delay due to the earth
                        st2t2s = stat_2_gcrs + tau*bet2s;
                        dGpe = 2*gme/c^3*...
                            log((norm(stat_1_gcrs)+rqu*stat_1_gcrs)/(norm(st2t2s)+rqu*st2t2s));
                        
                        tau = tau + dGps + dGpe;
                        
                        % aberrated source vector (yearly aberration)
                        s0d  = rqu;
                        coas = 1/(gam*(1+bet'*s0d'));
                        coab = coas*((gam-1)*(bet'*s0d')/(bet'*bet)+gam);
                        k1a  = coas*s0d'+coab*bet;
                        k2a  = k1a;
                        
                    case 3 % SSB acceleration
                        % calculation of the time delay following the consensus model
                        % + Titov 2011
                        % ------------------------------------------------------------

                        %(1) barycentric station vector (eq. 6)
                        xb1 = earth + stat_1_gcrs;
                        xb2 = earth + stat_2_gcrs;

                        %(2)&(3) differential gravitational delay due to celestial bodies
                        [Tgrav, pGammaSun] = grav_delay(xb1,xb2,vearth,b_gcrs,rqu,ephem,isc,opt);

                        %(4) differential gravitational delay due to the earth
                        Tgrave = 2*gme/c^3*...
                            log((norm(stat_1_gcrs)+rqu*stat_1_gcrs)/(norm(stat_2_gcrs)+rqu*stat_2_gcrs)); % (eq. 1)

                        %(5) total differential gravitational delay (eq. 2)
                        Tgrav = Tgrave + Tgrav;

                        %(6) vacuum delay
                        U     = gms/norm(sun); % gravitational potential at the geocenter
                        gamma = 1;
                        fac1  = (rqu*b_gcrs)/c;
                        term1 = 1-(1+gamma)*U/c^2-(norm(vearth))^2/(2*c^2)-(vearth'*v2)/c^2;
                        fac2  = ((vearth+acc*delt_accSSB)'*b_gcrs)/c^2;
                        term2 = 1+(rqu*vearth)/(2*c);
                        quot  = 1+(rqu*(vearth+v2+acc*delt_accSSB))/c;
                        tau   = (Tgrav - fac1*term1 - fac2*term2)/quot;          %(eq. 9)

                        pGammaSun = pGammaSun/quot;

                        %(7) aberrated source vector (eq. 15)
                        k1a = rqu + (vearth+v1)'/c - rqu*((rqu*(vearth+v1))')/c;
                        k2a = rqu + (vearth+v2)'/c - rqu*((rqu*(vearth+v2))')/c;
                        
                end %case
                
                % ##### Spacecrafts #####
            case 's'
                
                
                
                switch(del_mod_s)
                    
                    case 1 % Light time equation
                        
                        % Init.:
                        pGammaSun = []; % Not yet calculated for satelliet scans

                        if (stat_1_id == 1) || ~flag_fix_sat_pos_to_stat1                 
                            % Get reference epoch for interpolation of SC pos. + vel.:
                            ref_ind     = find((sources.s(scan(isc).iso).day == scan(isc).tim(3)) & (sources.s(scan(isc).iso).sec_of_day > sec_of_day));
                            if ref_ind(1) < 8
                                error('S.C. ephemeris data do not cover the required time (earlier epochs needed)! Add missing data to S.C. ephem. file!');
                            elseif ref_ind(1)+6 > length(sources.s(scan(isc).iso).day)
                                error('S.C. ephemeris data do not cover the required time (later epochs needed)! Add missing data to S.C. ephem. file!');
                            end
                            ref_ind             = ref_ind(1)-7 : 1 : ref_ind(1)+6; % suitable for interpolation with "lagint9.m" and "dt" < 1 sec
                            t_ref_sec_interpol  = sources.s(scan(isc).iso).sec_of_day(ref_ind);
                            t_integer_mjd       = floor(sources.s(scan(isc).iso).mjd(ref_ind));
                            t_ref_mjd           = t_integer_mjd(1);
                            offset_sec          = (t_integer_mjd - t_ref_mjd) * 86400; % full days since first interpolation epoch in [sec]
                            t_ref_sec_interpol  = t_ref_sec_interpol + offset_sec; % Add since first interpolation epoch in [sec]
                            t_ref_sec           = t_ref_sec_interpol(1);
                            t_ref_sec_interpol  = t_ref_sec_interpol - t_ref_sec;

                            t_integer_mjd_obs   = floor(mjd);
                            t_ref_offset_obs    = (t_integer_mjd_obs - t_ref_mjd) * 86400;
                            t_ref_sec_obs       = sec_of_day + t_ref_offset_obs;
                            t_ref_sec_obs       = t_ref_sec_obs -  t_ref_sec;



                            % Get spacecraft position at time of observation (CRF):
                            sc_x_crf = lagint9(t_ref_sec_interpol, sources.s(scan(isc).iso).x_crf(ref_ind), t_ref_sec_obs);
                            sc_y_crf = lagint9(t_ref_sec_interpol, sources.s(scan(isc).iso).y_crf(ref_ind), t_ref_sec_obs);
                            sc_z_crf = lagint9(t_ref_sec_interpol, sources.s(scan(isc).iso).z_crf(ref_ind), t_ref_sec_obs);
                            sc_pos_crf = [sc_x_crf; sc_y_crf; sc_z_crf];                                                        % (3,1), [m]


                            % Get spacecraft velocity at time of observation (CRF):
                            % - Calculated from positions one sc before and after the current obs. epoch
                            %    => Cont. for 2 sec => no further interpolation done
                            if ~sources.s(scan(isc).iso).flag_v_trf
                                sc_x_crf_m  = lagint9(t_ref_sec_interpol, sources.s(scan(isc).iso).x_crf(ref_ind), t_ref_sec_obs - 1);
                                sc_y_crf_m  = lagint9(t_ref_sec_interpol, sources.s(scan(isc).iso).y_crf(ref_ind), t_ref_sec_obs - 1);
                                sc_z_crf_m  = lagint9(t_ref_sec_interpol, sources.s(scan(isc).iso).z_crf(ref_ind), t_ref_sec_obs - 1);
                                sc_x_crf_p  = lagint9(t_ref_sec_interpol, sources.s(scan(isc).iso).x_crf(ref_ind), t_ref_sec_obs + 1);
                                sc_y_crf_p  = lagint9(t_ref_sec_interpol, sources.s(scan(isc).iso).y_crf(ref_ind), t_ref_sec_obs + 1);
                                sc_z_crf_p  = lagint9(t_ref_sec_interpol, sources.s(scan(isc).iso).z_crf(ref_ind), t_ref_sec_obs + 1);
                                sc_vel_crf  = ([sc_x_crf_p; sc_y_crf_p; sc_z_crf_p] - [sc_x_crf_m; sc_y_crf_m; sc_z_crf_m])/2;      %(3,1), [m/sec]
                            end



                            % Get spacecraft position at the time of emission (CRF):
                            % - According to approach "geocneu = 0" in vie_mod_tie.m (lines 735-763) by L. Plank
                            sc_pos_crf_tmp(1, :) = sc_pos_crf';
                            t_ref_sec_obs_tmp = t_ref_sec_obs;

                            % Iteration init.:
                            number_of_iterations    = 0;
                            dt                      = 999999;
                            ddt                     = 999999;

                            while(abs(ddt) > ddt_threshold) 
                                number_of_iterations = number_of_iterations + 1;

                                % Get spacecraft velocity at time of observation (CRF):
                                % - from ephem. file, if available
                                if sources.s(scan(isc).iso).flag_v_trf
                                    sc_vx_crf = lagint9(t_ref_sec_interpol, sources.s(scan(isc).iso).vx_crf(ref_ind), t_ref_sec_obs_tmp);
                                    sc_vy_crf = lagint9(t_ref_sec_interpol, sources.s(scan(isc).iso).vy_crf(ref_ind), t_ref_sec_obs_tmp);
                                    sc_vz_crf = lagint9(t_ref_sec_interpol, sources.s(scan(isc).iso).vz_crf(ref_ind), t_ref_sec_obs_tmp);
                                    sc_vel_crf = [sc_vx_crf; sc_vy_crf; sc_vz_crf];
                                end

                                % Correction:
                                dt_old  = dt;
                                dt      = norm(sc_pos_crf - stat_1_gcrs)/c - ((sc_pos_crf' - stat_1_gcrs')*sc_vel_crf)/c^2; % [sec] (6.4)
                                t_ref_sec_obs_tmp = t_ref_sec_obs - dt; % correctd epoch [sec] since "ref. time" (t_ref_sec, t_ref_mjd)

                                %iteration1
                                sc_x_crf = lagint9(t_ref_sec_interpol, sources.s(scan(isc).iso).x_crf(ref_ind), t_ref_sec_obs_tmp);
                                sc_y_crf = lagint9(t_ref_sec_interpol, sources.s(scan(isc).iso).y_crf(ref_ind), t_ref_sec_obs_tmp);
                                sc_z_crf = lagint9(t_ref_sec_interpol, sources.s(scan(isc).iso).z_crf(ref_ind), t_ref_sec_obs_tmp);
                                sc_pos_crf = [sc_x_crf; sc_y_crf; sc_z_crf];
                                sc_pos_crf_tmp(number_of_iterations +1, :) = sc_pos_crf';

                                ddt = dt_old - dt;
                                if number_of_iterations >= max_iterations
                                    fprintf(' Warning: Max. number of iterations (%d) for near filed delay reached! ddt = %f sec', max_iterations, ddt);
                                    break;
                                end

                            end
                        end

                        
                        % Vector station-spacecraft (source vectors) at the time of signal emission (spacecraft) and reception at station one (stations)
                        L1  = sc_pos_crf' - stat_1_gcrs';  % (1x3)
                        L2  = sc_pos_crf' - stat_2_gcrs';
                                
                        % gravitational potential @geocentre / all except earth
                        sunb = sun + earth;
                        Wsun  = ephem.gms/norm(earth-sunb);
                        moonb = moon + earth;
                        Wmoon = gmm/norm(earth-moonb);
                        Wmerc = ephem.gmmerc/norm(earth-ephem.merc(isc).xbar);
                        Wvenu = ephem.gmvenu/norm(earth-ephem.venu(isc).xbar);
                        Wmars = ephem.gmmars/norm(earth-ephem.mars(isc).xbar);
                        Wjupi = ephem.gmjupi/norm(earth-ephem.jupi(isc).xbar);
                        Wsatu = ephem.gmsatu/norm(earth-ephem.satu(isc).xbar);
                        Wuran = ephem.gmuran/norm(earth-ephem.uran(isc).xbar);
                        Wnept = ephem.gmnept/norm(earth-ephem.nept(isc).xbar);
                        Wplut = ephem.gmplut/norm(earth-ephem.plut(isc).xbar);
                        We = Wplut + Wnept + Wuran + Wsatu + Wjupi + Wmars + Wvenu + Wmerc + Wmoon + Wsun;
                                
                        du0     = (norm(L2) - norm(L1))/c; % Difference in travel time not considering retarded BL effect, [sec]
                        n       = L2/norm(L2);
                        dugr    = 2*gme/c^3*log(...
                                    ((norm(stat_2_gcrs) + norm(sc_pos_crf) + norm(stat_2_gcrs - sc_pos_crf)) * (norm(stat_1_gcrs) + norm(sc_pos_crf) - norm(stat_1_gcrs - sc_pos_crf)))...
                                    /((norm(stat_2_gcrs) + norm(sc_pos_crf) - norm(stat_2_gcrs - sc_pos_crf)) * (norm(stat_1_gcrs) + norm(sc_pos_crf) + norm(stat_1_gcrs - sc_pos_crf)))); % Gravitaional effect on travel time [sec]
                            
                        % For the GEOCENTR dugr becomes -inf!
                        if strcmp(antenna(stat_1_id).name,'GEOCENTR') || strcmp(antenna(stat_2_id).name,'GEOCENTR')
                            dugr = 0;
                            du  = du0*(1-n*v2/c);% + dugr - 1/c^2*((v2'*v2)/2+We) * du0; % Delay:  Klioner, 1991, formular (6.3) [sec]
                        else
                            du  = du0*(1-n*v2/c) + dugr - 1/c^2*((v2'*v2)/2+We) * du0; % Delay:  Klioner, 1991, formular (6.3) [sec]
                        end
                                
                        K   = (L1+L2)/(norm(L1)+norm(L2)); % mittlere Richtungsvektor
%                         b   = stat_1_gcrs - stat_2_gcrs; % baseline vector
%                         phi = acos(K * b /(norm(K) * norm(b))); %[rad] % 
                                
%                         dphi    = -(norm(K) * norm(b)) * sin(phi); % 
%                         tauobs  = norm(K) * norm(b) / c*cos(phi);
%                         delphi  = phi;

                        % Vacuum delay:
                        tau = du;

                        % Source vectors:
                        k1a = L1; % Source vector stat 1 
                        k2a = L2; % Source vector stat 2 
                        
                        rqu = K;
                        rq  = K;
                        
                        % #### partial derivative of the delay w.r.t.:  ####
                        % Position of the space craft:

                        % Partial derivative of du instead of du0 (without gravitational potential, with retarded BL correction)
                        % dtdsp1 = -(L1)/norm(L1)+L2/norm(L2); % Ableiutng von du0 nach ws
                        % dtdsp2 = (v2'/c) .* ( (L2/norm(L2) - L1/norm(L1)) .* L2/norm(L2)   +   ((norm(L2) + sc_pos_crf' .* (L2 / (2*norm(L2)^4))) / norm(L2)^2) .* (norm(L2) - norm(L1)) );
                        % dtds0 = dtdsp1-dtdsp2; % In GCRF
                        
                        
                        % Neu:
                        nL1 = norm(L1);
                        nL2 = norm(L2);
                                                
                        % PD in GCRF:
                        dudws1_part1 = (sc_pos_crf(1) - stat_2_gcrs(1)) / nL2    -   (sc_pos_crf(1) - stat_1_gcrs(1)) / nL1; 
                        dudws2_part1 = (sc_pos_crf(2) - stat_2_gcrs(2)) / nL2    -   (sc_pos_crf(2) - stat_1_gcrs(2)) / nL1; 
                        dudws3_part1 = (sc_pos_crf(3) - stat_2_gcrs(3)) / nL2    -   (sc_pos_crf(3) - stat_1_gcrs(3)) / nL1; 
                        
                        dudws1_part2 = v2(1)   -   (sc_pos_crf(1) - stat_1_gcrs(1))/nL1 * 1/nL2 * (L2*v2)    +    nL1 * (sc_pos_crf(1) - stat_2_gcrs(1)) * 1/ nL2^3 * (L2*v2)   -  nL1 * 1/nL2 * v2(1); 
                        dudws2_part2 = v2(2)   -   (sc_pos_crf(2) - stat_1_gcrs(2))/nL1 * 1/nL2 * (L2*v2)    +    nL1 * (sc_pos_crf(2) - stat_2_gcrs(2)) * 1/ nL2^3 * (L2*v2)   -  nL1 * 1/nL2 * v2(2); 
                        dudws3_part2 = v2(3)   -   (sc_pos_crf(3) - stat_1_gcrs(3))/nL1 * 1/nL2 * (L2*v2)    +    nL1 * (sc_pos_crf(3) - stat_2_gcrs(3)) * 1/ nL2^3 * (L2*v2)   -  nL1 * 1/nL2 * v2(3); 
                        
                        dudws1 = dudws1_part1 / c   -   dudws1_part2 / c^2; % [sec/m]
                        dudws2 = dudws2_part1 / c   -   dudws2_part2 / c^2; % [sec/m]
                        dudws3 = dudws3_part1 / c   -   dudws3_part2 / c^2; % [sec/m]
                        
                        pd_satpos_gcrf = [dudws1; dudws2; dudws3];          % PD of delay time du w.r.t. satellite pos. in GCRF [sec/m], (3x1 vector) 
                        pd_satpos_gcrf = pd_satpos_gcrf * c; % * 100 / 100; % Unit conversion: [sec/m] => [cm/cm] = []; => estimates will be in [cm]
                        
                        % Rotation of PD to RSW system:
                        [~, ~, transmat] = rv2rsw(sc_pos_crf, sc_vel_crf);
                        pd_satpos_rsw = transmat*pd_satpos_gcrf;
                                                
                        % Rotate PD vector to TRF and RSW system:
                        % .... add code here!
                        
                        % +++++ Test PD w.r.t. sat.Pos +++++
                        % Nummerical derivation of du0 w.r.t. ws1 (du0/dws1)
                        flag_calc_numerical_derivation_of_sat_pos = false;
                        if flag_calc_numerical_derivation_of_sat_pos
                            sc_pos_crf_plus = sc_pos_crf + [0.5; 0; 0]; % plus 0.5 m in ws1
                            L1_plus     = sc_pos_crf_plus' - stat_1_gcrs';
                            L2_plus     = sc_pos_crf_plus' - stat_2_gcrs';
                            sc_pos_crf_minus = sc_pos_crf - [0.5; 0; 0]; % minus 0.5 m in ws1
                            L1_minus    = sc_pos_crf_minus' - stat_1_gcrs';
                            L2_minus     = sc_pos_crf_minus' - stat_2_gcrs';
                            dudws1_part1_num     = ( ((norm(L2_plus) - norm(L1_plus))) - ((norm(L2_minus) - norm(L1_minus))) ) / 1;
                            n_plus      =  L2_plus / norm(L2_plus);
                            n_minus     =  L2_minus / norm(L2_minus);
                            dudws1_part2_num = ((norm(L2_plus) - norm(L1_plus)) * (n_plus * v2))   -   ((norm(L2_minus) - norm(L1_minus)) * (n_minus * v2));               

                            sc_pos_crf_tmp1 = sc_pos_crf + [0; 0.5; 0]; % plus 0.5 m in ws1
                            L1_plus     = sc_pos_crf_tmp1' - stat_1_gcrs';
                            L2_plus     = sc_pos_crf_tmp1' - stat_2_gcrs';
                            sc_pos_crf_tmp1 = sc_pos_crf - [0; 0.5; 0]; % minus 0.5 m in ws1
                            L1_minus    = sc_pos_crf_tmp1' - stat_1_gcrs';
                            L2_minus     = sc_pos_crf_tmp1' - stat_2_gcrs';
                            dudws2_part1_num     = ( ((norm(L2_plus) - norm(L1_plus))) - ((norm(L2_minus) - norm(L1_minus))) ) / 1;
                            n_plus      =  L2_plus / norm(L2_plus);
                            n_minus     =  L2_minus / norm(L2_minus);
                            dudws2_part2_num = ((norm(L2_plus) - norm(L1_plus)) * (n_plus * v2))   -   ((norm(L2_minus) - norm(L1_minus)) * (n_minus * v2));

                            sc_pos_crf_tmp1 = sc_pos_crf + [0; 0; 0.5]; % plus 0.5 m in ws1
                            L1_plus     = sc_pos_crf_tmp1' - stat_1_gcrs';
                            L2_plus     = sc_pos_crf_tmp1' - stat_2_gcrs';
                            sc_pos_crf_tmp1 = sc_pos_crf - [0; 0; 0.5]; % minus 0.5 m in ws1
                            L1_minus    = sc_pos_crf_tmp1' - stat_1_gcrs';
                            L2_minus     = sc_pos_crf_tmp1' - stat_2_gcrs';
                            dudws3_part1_num     = ( ((norm(L2_plus) - norm(L1_plus))) - ((norm(L2_minus) - norm(L1_minus))) ) / 1;
                            n_plus      =  L2_plus / norm(L2_plus);
                            n_minus     =  L2_minus / norm(L2_minus);
                            dudws3_part2_num = ((norm(L2_plus) - norm(L1_plus)) * (n_plus * v2))   -   ((norm(L2_minus) - norm(L1_minus)) * (n_minus * v2));
                        end
                        % ==> analytical and numerical solutions are equal!
                        % ----- Test PD w.r.t. sat.Pos -----          
                                                
                        % Rotation to sat.-ref.sys.
                        % =>> Seems to be wrong! See: Vallado, Fundamanetzals of Astrodynamics and Applications, p.163, RSW system!!!
                        % crosso = cross(sc_pos_crf',sc_vel_crf);                     % Get cross-track direction (unit vector)
                        % dtds(1)=dtds0(1:3)*(sc_pos_crf'/norm(sc_pos_crf'))';    % radial
                        % dtds(2)=dtds0(1:3)*(sc_vel_crf/norm(sc_vel_crf));               % along track (falsch! => Hier geh򲴺 cross(cross-track, radial), NICHT genau in v-Richtung!
                        % dtds(3)=dtds0(1:3)*(crosso/norm(crosso))';              % cross track
                        
                        % Station coordinates:
                        % - P.D. of du0 w.r.t. station coordinates (in TRF!)
                        ps1 = -t2c'*(L1'/norm(L1)); % Partial derivative: station 1
                        ps2 = t2c'*(L2'/norm(L2));  % Partial derivative: station 2
                        
                end % switch(del_mod_s)
                
        end % switch(scan(isc).obs_type)
        
        
        % the following corrections are the same for both models (Sekido & Fukushima, p.141)
        % ---- delay calculation finished -------
        % ---- further corrections: -------
        for kstat=1:2
            switch kstat
                case 1
                    stnum = stat_1_id; ka = k1a; % stt = stt1;
                case 2
                    stnum = stat_2_id; ka = k2a; % stt = stt2;
            end
            
            % only do once for one station for one scan
            if isempty(scan(isc).stat(stnum).zd)
                
                % Get corrected station coordinates for actual observation epoch:
                [phi, lam, hell] = xyz2ell(scan(isc).stat(stnum).x);
                
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
                
                
                if strcmpi(parameter.vie_init.tropSource.name,'indModeling')   ||   isempty(scan(isc).stat(stnum).trop)

                    % TROPOSPHERE
                    
                    % apriori zenith delay Lz [m] Marini tropospheric model

                    % hydrostatic
                    if strcmpi(parameter.vie_init.zhd,'in situ')
                        pres = scan(isc).stat(stnum).pres;  % [hPa]
                        zdry = 0.0022768*pres / (1-0.00266*cos(2*phi)-(0.28e-6*hell));   %[m]
                    elseif strcmpi(parameter.vie_init.zhd,'no')
                        zdry = 0;
                    elseif strcmpi(parameter.vie_init.zhd,'vmf3')
                        zdry = scan(isc).stat(stnum).zhdt;
                    elseif strcmpi(parameter.vie_init.zhd,'vmf1')
                        zdry = scan(isc).stat(stnum).zhdt;
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
                        e      = scan(isc).stat(stnum).e;  % [hPa]
                        Tm     = antenna(stnum).gpt3.Tm;
                        lambda = antenna(stnum).gpt3.lambda;
                        zwet   = asknewet ( e , Tm , lambda );
                    elseif strcmpi(parameter.vie_init.zwd,'vmf3')
                        zwet   = scan(isc).stat(stnum).zwdt;
                    elseif strcmpi(parameter.vie_init.zwd,'vmf1')
                        zwet   = scan(isc).stat(stnum).zwdt;
                    elseif strcmpi(parameter.vie_init.zwd,'gpt3')
                        e      = antenna(stnum).gpt3.e;
                        Tm     = antenna(stnum).gpt3.Tm;
                        lambda = antenna(stnum).gpt3.lambda;
                        zwet   = asknewet ( e , Tm , lambda );
                    else
                        error('Something is wrong here...');
                    end
                    
                    
                    % mapping function
                    
                    % hydrostatic
                    aht = scan(isc).stat(stnum).aht;
                    if strcmpi(parameter.vie_mod.mfh,'vmf3') == 1 && aht ~= 0
                        [mfh,~] = vmf3(aht,awt,mjd,phi,lam,zd);   % VMF3
                    elseif strcmpi(parameter.vie_mod.mfh,'vmf1') == 1 && aht ~= 0
                        [mfh,~] = vmf1(aht,awt,mjd,phi,zd);   % VMF1
                    elseif strcmpi(parameter.vie_mod.mfh,'gpt3')
                        [~,~,~,~,~,aht,awt,~,~,~,~,~,~] = gpt3_5_fast (mjd,phi,lam,hell,0,cell_grid_GPT3);
                        [mfh,~] = vmf3_ht (aht,awt,mjd,phi,lam,hell,zd);   % GPT3
                    else
                        [mfh,~] = vmf3_ht (antenna(stnum).gpt3.ah,antenna(stnum).gpt3.aw,mjd,phi,lam,hell,zd);   % GPT3 as backup if aht = 0
                    end
                    
                    % wet
                    awt = scan(isc).stat(stnum).awt;
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
                    scan(isc).stat(stnum).mfw   = mfw;   % used in vie_lsm as partial
                    scan(isc).stat(stnum).zdry  = zdry;
                    scan(isc).stat(stnum).zwet  = zwet;
                    scan(isc).stat(stnum).aoalt   = aoalt; % antenna axis offset altitude correction [m]                    
                    scan(isc).stat(stnum).trop  = ((zdry+aoalt) * mfh + zwet * mfw + aprgrd)/c; % contains the whole (asymmetric) delay [sec]
 
                end 
                
                
                % THERMAL DEFORMATION (Haas et al.,1998; Skurihina 2000)
                if parameter.vie_mod.therm == 1
                    thermal = antenna(stnum).thermal;
                    temp = scan(isc).stat(stnum).temp;  % measured temperature [C] or from gpt3
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
                        temp = scan(isc).stat(stnum).temp;
                        delay_temp = (temp - 19) * 0.47 * 1e-3;  % [m]
                        delay_temp = delay_temp / c;  % [m] --> [s]
                        
                        gravdef_corr = gravdef_corr + delay_temp;
                    end                    
                else
                    gravdef_corr = 0;
                end
                
                % store in scan
                scan(isc).stat(stnum).axkt  = axkt;
                scan(isc).stat(stnum).therm = therm_d;
                scan(isc).stat(stnum).az    = azim;
                scan(isc).stat(stnum).zd    = zd; % corrected for aberration, not for refraction!
                scan(isc).stat(stnum).daxkt  = daxkt; % [-] partial derivative AO
                scan(isc).stat(stnum).gravdef = gravdef_corr;

                
            end % if zd is empty
        end % kstat 1:2
        
        %(8) add geometric part of tropospheric propagation delay
        atm1  = scan(isc).stat(stat_1_id).trop;  %[sec]
        atm2  = scan(isc).stat(stat_2_id).trop;  %[sec]
        
        % For the GEOCENTR atm1 becomes -inf!
        % => zenith delay form second station is taken instead as a prelimiary solution
        if strcmp(antenna(stat_1_id).name,'GEOCENTR')
            atm1 = scan(isc).stat(stat_2_id).zdry/c;
            scan(isc).stat(stat_1_id).trop = atm1;
        elseif strcmp(antenna(stat_2_id).name,'GEOCENTR')
            atm2 = scan(isc).stat(stat_1_id).zdry/c;
            scan(isc).stat(stat_2_id).trop = atm2;
        end
        
        tpd_g = atm1*(rqu*(v2-v1))/c; % 
        tau   = tau + tpd_g;         %(eq. 11)
        
        %(9) total delay
        c_trop = atm2 - atm1;
        tau    = tau + c_trop; %(eq. 12)
        
        % further corrections:
        % axis offset correction
        c_axis  = scan(isc).stat(stat_2_id).axkt - scan(isc).stat(stat_1_id).axkt; %[sec]
        
        % thermal deformation (Attention: Station1 - Station2)
        c_therm = scan(isc).stat(stat_1_id).therm - scan(isc).stat(stat_2_id).therm; % [sec]
        
        % gravitational deformation correction
        c_gravdef = scan(isc).stat(stat_2_id).gravdef - scan(isc).stat(stat_1_id).gravdef; % [sec]
        
        % add
        tau = c_axis + c_therm + c_gravdef + tau; % [sec]
        
        % + EXTERNAL IONOSPERIC DELAY +
        if strcmp(parameter.vie_init.iono, 'ext')
            % use iono delay from external file and add it
            
            % only do the first time
            if ~exist('iondata', 'var')
                fprintf('Start loading external ionospheric file\n');
                [iondata, ionFileFoundLog] = load_ionfile(parameter,session);
            end
            
            % if iono file was found -> apply correction
            if ionFileFoundLog == 1
                
                % find two stationnames of current observations
                statNameI1      = antenna(scan(isc).obs(iobs).i1).name;
                statNameI2      = antenna(scan(isc).obs(iobs).i2).name;
                [dion1, dion2]  = get_iondel(iondata,scan(isc).tim,statNameI1,statNameI2);
                ionoDel         = dion2 - dion1; %[sec]
                
                % add to observed delay
                scan(isc).stat(stat_1_id).iono  = dion1;
                scan(isc).stat(stat_2_id).iono  = dion2;
                scan(isc).obs(iobs).ionDelext   = ionoDel;
                scan(isc).obs(iobs).obs         = scan(isc).obs(iobs).obs - ionoDel;
            else
                warning('Ion. correction not found in external file! Correction set to zero.\n');
                scan(isc).stat(stat_1_id).iono  = 0;
                scan(isc).stat(stat_2_id).iono  = 0;
                scan(isc).obs(iobs).ionDelext   = 0;
                
            end
        end
        %  - EXTERNAL IONOSPERIC DELAY -
        
        % SOURCE STRUCTURE +
        if parameter.vie_mod.ssou==1 || parameter.vie_mod.write_jet==1
            stat_1_gcrs = scan(isc).obs(iobs).i1;
            stat_2_gcrs = scan(isc).obs(iobs).i2;
            ind=strmatch(sources.q(scan(isc).iso).name,cat_comp.name,'exact');
            if isempty(ind)
                nam2=sounamivs2iers(sources.q(scan(isc).iso).name);
                ind = strmatch(nam2,cat_comp.name,'exact');
            end
            if isempty(ind)
                disp(strcat('source ',sources.q(scan(isc).iso).name,' not found in ss catalogue.'));
                soucorr=0; jetang=90; jetjb=0; uvrange=0; uu=0; vv=0;
            else
                sources.q(scan(isc).iso).sou_model=[cat_comp.flux(ind),cat_comp.maj(ind),cat_comp.min(ind),cat_comp.angle(ind),cat_comp.dRA(ind),cat_comp.dDec(ind)];
                [soucorr,uu,vv]=modDelay(sources.q(scan(isc).iso).sou_model,scan(isc).stat(stat_1_gcrs).x,scan(isc).stat(stat_2_gcrs).x,([8213 8252 8353 8513 8733 8853 8913 8933]+4), 8217, sources.q(scan(isc).iso).ra2000, sources.q(scan(isc).iso).de2000, scan(isc).mjd);
                soucorr=soucorr*1e-12;
                jetvec=(sources.q(scan(isc).iso).sou_model(2,5:6));
                jetvec=jetvec/norm(jetvec);
                uvvec=[uu;vv];
                uvvec=uvvec/norm(uvvec);
                jetang=acos(abs(jetvec*uvvec))*180/pi;
                uvrange=(jetvec)*[uu;vv];
                jetjb=(90-jetang)*(norm(b_gcrs)/6371000);
            end
        end
        
        % correct for source structure
        if parameter.vie_mod.ssou==1
            tau=tau+soucorr;
        end
        
        % writing jetang to external file
        if parameter.vie_mod.write_jet==1
            fprintf(fidjet,'%s %s %4.12f %4.12f\n',antenna(scan(isc).obs(iobs).i1).name,antenna(scan(isc).obs(iobs).i2).name,mjd,jetang);
            fprintf(fidjetuv,'%s %s %4.12f %4.12f\n',antenna(scan(isc).obs(iobs).i1).name,antenna(scan(isc).obs(iobs).i2).name,mjd,uvrange);
            fprintf(fidjetjb,'%s %s %4.12f %4.12f\n',antenna(scan(isc).obs(iobs).i1).name,antenna(scan(isc).obs(iobs).i2).name,mjd,jetjb);
        end
        % SOURCE STRUCTURE -
        
        
        % ������������������������������������������������������������������������
        % �  4. PARTIAL DERIVATIVES                                              �
        % ������������������������������������������������������������������������
        %
        % according to MODEST Handbook 1994 p.46 ff
        beta = vearth/c;
        b2   = (v2 + vearth)/c;
        gam  = 1/sqrt(1-beta'*beta);
        rho  = 1+rq*b2;
        dij  = eye(3);
        
        psi = -(gam*(1-beta'*b2)*rq'/rho+gam*beta);               %[3,1]  (2.220)
        E   = dij +((gam-1)*beta/(beta'*beta)-gam*b2)*beta';      %[3,3]  (2.243)
        K   = E*psi;                                              %[3,1]  (2.264)
        B   = K'*t2c;                                             %[1,3]  (2.249 with t2c = Qjk)
        M   = (dij -(rq'*b2')/rho)*(-gam*(1-b2'*beta)*(E*b_gcrs)/rho); %[3,1]
        
        % wrt EOP per baseline [m]
        pdx = K'* (dQdxp * b_trs')/c;
        pdy = K'* (dQdyp * b_trs')/c;
        put = K'* (dQdut * b_trs')/c;
        pdX = K'* (dQddX * b_trs')/c;
        pdY = K'* (dQddY * b_trs')/c;
        
        % wrt source coordinates [cm/mas]
        % RaDec position of quasars
        % - set =0, if no quasar was observed in this scan (e.g. scan to spacecraft)
        if strcmp(scan(isc).obs_type, 'q')
            parra = (drqdra * M)*pi()/180/3600000*100;  % (2.230)
            parde = (drqdde * M)*pi()/180/3600000*100;  % (2.231)  % [cm/mas]
			psou = [parra, parde];
        else 
%             parra = 0;
%             parde = 0;  % [cm/mas]
            psou = [];
        end
        
        % wrt station coordinates
        if strcmp(scan(isc).obs_type, 'q')
            ps1  = -B; % (2.248)
            ps2  =  B; % (2.248)
        end
        
        % Axis Offset
        paxkt_st1 = -scan(isc).stat(stat_1_id).daxkt; %[-]
        paxkt_st2 =  scan(isc).stat(stat_2_id).daxkt; %[-]
        
        
        % Amplitudes of the seasonal variation of station position
        pAcr_st1 = (- scan(isc).stat(stat_1_id).pAcr_xyz)* B';
        pAce_st1 = (- scan(isc).stat(stat_1_id).pAce_xyz)* B';
        pAcn_st1 = (- scan(isc).stat(stat_1_id).pAcn_xyz)* B';
        pAsr_st1 = (- scan(isc).stat(stat_1_id).pAsr_xyz)* B';
        pAse_st1 = (- scan(isc).stat(stat_1_id).pAse_xyz)* B';
        pAsn_st1 = (- scan(isc).stat(stat_1_id).pAsn_xyz)* B';
        
        pAcr_st2 = (scan(isc).stat(stat_2_id).pAcr_xyz)* B';
        pAce_st2 = (scan(isc).stat(stat_2_id).pAce_xyz)* B';
        pAcn_st2 = (scan(isc).stat(stat_2_id).pAcn_xyz)* B';
        pAsr_st2 = (scan(isc).stat(stat_2_id).pAsr_xyz)* B';
        pAse_st2 = (scan(isc).stat(stat_2_id).pAse_xyz)* B';
        pAsn_st2 = (scan(isc).stat(stat_2_id).pAsn_xyz)* B';
        
        % pole tide Love & Shida  per baseline
        phpole_bl = (scan(isc).stat(stat_2_id).phpole - scan(isc).stat(stat_1_id).phpole)* B';
        plpole_bl = (scan(isc).stat(stat_2_id).plpole - scan(isc).stat(stat_1_id).plpole)* B';
        
        
        % APL RG
        prg_st1 = (- scan(isc).stat(stat_1_id).prg)* B';
        prg_st2 = (scan(isc).stat(stat_2_id).prg)* B';
        
        
        % tidal ERP partials
        % wrt tidal ERP term per baseline [m]
        
        if opt.est_erpsp == 1
            pap = zeros(1,length(cosag));
            pam = zeros(1,length(cosag2));
            pbp = zeros(1,length(cosag));
            pbm = zeros(1,length(cosag2));
            puc = zeros(1,length(cosag));
            pus = zeros(1,length(cosag));
            
            % prograde polar motion and dUT1 terms
            for tidn = 1 : length(cosag)
                pap(tidn) = K'* ((PNn*Rall* (R3ss*(dR2xp.*(-cosag(tidn))*R1yp + ...
                    R2xp*dR1yp.*sinag(tidn)))) * b_trs')/c;
                pbp(tidn) = K'* ((PNn*Rall* (R3ss*(dR2xp.*sinag(tidn)*R1yp + ...
                    R2xp*dR1yp.*cosag(tidn)))) * b_trs')/c;
                puc(tidn) = K'* ((PNn    * dRdutn.*cosag(tidn) * Wn) * b_trs')/c;
                pus(tidn) = K'* ((PNn    * dRdutn.*sinag(tidn) * Wn) * b_trs')/c;
            end
            
            % retrograde semidiurnal polar motion
            for tidn = 1 : length(cosag2)
                pam(tidn) = K'* ((PNn*Rall* (R3ss*(dR2xp.*(-cosag2(tidn))*R1yp + ...
                    R2xp*dR1yp.*(-sinag2(tidn))))) * b_trs')/c;
                pbm(tidn) = K'* ((PNn*Rall* (R3ss*(dR2xp.*(-sinag2(tidn))*R1yp + ...
                    R2xp*dR1yp.*cosag2(tidn)))) * b_trs')/c;
            end
        end
        
        if opt.est_FCNnut == 1
            pFCNnut = K'* (dQdfcn * b_trs')/c;
        end
        if opt.est_FCNnutAmpl == 1
            pFCNnutAc = K'* (dQdfcnAc * b_trs')/c; %[s]
            pFCNnutAs = K'* (dQdfcnAs * b_trs')/c; %[s]
        end
        
        % Love & Shida per baseline
        pLove_bl  = (scan(isc).stat(stat_2_id).pLove  - scan(isc).stat(stat_1_id).pLove)* B';
        pShida_bl = (scan(isc).stat(stat_2_id).pShida - scan(isc).stat(stat_1_id).pShida)* B';
        pFCN_bl   = (scan(isc).stat(stat_2_id).pFCN - scan(isc).stat(stat_1_id).pFCN)* B';
        
        % wrt SBB acceleration
        pacc = delt_accSSB/c^2*((rq*b_gcrs)*rq'-b_gcrs)...
            -delt_accSSB/c^3*((rq*(vearth+v2))*b_gcrs...
            +((rq*vearth)*b_gcrs)/2 - (b_gcrs'*vearth)*rq'); %[sec^3/m]
        pacc = pacc'./100; %[sec^3/cm]
        
        
        
        
        % ������������������������������������������������������������������������
        % �  5. STORE RESULTS                                                    �
        % ������������������������������������������������������������������������
        
        scan(isc).obs(iobs).com = tau; %[sec]
        
        % partial derivatives of the delay  w.r.t.: 
        scan(isc).obs(iobs).psou                = psou;             % source coordinates [cm/mas]; =[0,0], if no quasar was observed
        scan(isc).obs(iobs).psou_sat_gcrf       = pd_satpos_gcrf;   % satellite coordinates in GCRF [sec/m]
        scan(isc).obs(iobs).psou_sat_trf        = pd_satpos_trf;    % satellite coordinates in TRF [sec/m]
        scan(isc).obs(iobs).psou_sat_rsw        = pd_satpos_rsw;    % satellite coordinates in RSW system ("satellite coord. sys.") [sec/m]
        scan(isc).obs(iobs).pnut                = [pdX,pdY];        % dX, dY [sec/rad]
        scan(isc).obs(iobs).ppol                = [pdx,pdy,put];    % xpol, ypol, dut1 [sec/rad]
        scan(isc).obs(iobs).pstat1              = ps1;              % station1
        scan(isc).obs(iobs).pstat2              = ps2;              % station2
        scan(isc).obs(iobs).pAO_st1             = paxkt_st1;        % Axis offset at station 1 [-]
        scan(isc).obs(iobs).pAO_st2             = paxkt_st2;        % Axis offset at station 2 [-]
        
        if parameter.vie_mod.ssou==1 || parameter.vie_mod.write_jet==1
            scan(isc).obs(iobs).jetang  =        jetang;
            scan(isc).obs(iobs).uvrange =       uvrange;
            scan(isc).obs(iobs).jetjb   =         jetjb;
            scan(isc).obs(iobs).uu      =            uu;
            scan(isc).obs(iobs).vv      =            vv;
            scan(isc).obs(iobs).soucorr =       soucorr; % [sec]
        end
        
        
        
        scan(isc).obs(iobs).pAcr_st1 = pAcr_st1;  % Amplitudes of the seasonal variation in station pos.
        scan(isc).obs(iobs).pAce_st1 = pAce_st1;  % [-]
        scan(isc).obs(iobs).pAcn_st1 = pAcn_st1;
        scan(isc).obs(iobs).pAsr_st1 = pAsr_st1;
        scan(isc).obs(iobs).pAse_st1 = pAse_st1;
        scan(isc).obs(iobs).pAsn_st1 = pAsn_st1;
        
        scan(isc).obs(iobs).pAcr_st2 = pAcr_st2;  % Amplitudes of the seasonal variation in station pos.
        scan(isc).obs(iobs).pAce_st2 = pAce_st2;  % [-]
        scan(isc).obs(iobs).pAcn_st2 = pAcn_st2;
        scan(isc).obs(iobs).pAsr_st2 = pAsr_st2;
        scan(isc).obs(iobs).pAse_st2 = pAse_st2;
        scan(isc).obs(iobs).pAsn_st2 = pAsn_st2;
        
        scan(isc).obs(iobs).phpole   = phpole_bl; % [cm] for pole tide Love number
        scan(isc).obs(iobs).plpole   = plpole_bl; % [cm] for pole tide Shida number
        
        scan(isc).obs(iobs).prg_st1  = prg_st1; %[hPa]  APL RG
        scan(isc).obs(iobs).prg_st2  = prg_st2; %[hPa]	APL RG
        
        % tidal ERP
        if opt.est_erpsp == 1
            scan(isc).obs(iobs).ppap  = pap;  % prograde pm cos
            scan(isc).obs(iobs).ppam  = pam;  % retrograde pm cos
            scan(isc).obs(iobs).ppbp  = pbp;  % prograde pm sin
            scan(isc).obs(iobs).ppbm  = pbm;  % retrograde pm sin
            scan(isc).obs(iobs).putc  = puc;  % ut1 cos
            scan(isc).obs(iobs).puts  = pus;  % ut1 sin
        end
        
        if opt.est_FCNnut ==1
            scan(isc).obs(iobs).pFCNnut     =  pFCNnut;     % [sec*day]
        else
            scan(isc).obs(iobs).pFCNnut     =  0;           % [sec*day]
        end
        if opt.est_FCNnutAmpl==1
            scan(isc).obs(iobs).pFCNnutAc   =  pFCNnutAc;   % [sec]
            scan(isc).obs(iobs).pFCNnutAs   =  pFCNnutAs;   % [sec]
        else
            scan(isc).obs(iobs).pFCNnutAc   =  0;           % [sec]
            scan(isc).obs(iobs).pFCNnutAs   =  0;           % [sec]
        end
        
        scan(isc).obs(iobs).pLove    = pLove_bl;    % [cm] for Love numbers
        scan(isc).obs(iobs).pShida   = pShida_bl;   % [cm] for Shida
        scan(isc).obs(iobs).pFCN     = pFCN_bl;     % [cm] for FCN
        scan(isc).obs(iobs).pacc     = pacc;        % dt/dacc SSB acceleration [sec^3/cm]
        scan(isc).obs(iobs).pGamma   = pGammaSun;   % [sec]       
        scan(isc).obs(iobs).pscale   = -fac1;       % [sec] Correction to the scale factor

    if flag_save_results
        i_result = i_result + 1;
        result.com(i_result)   = tau; % Computed delay
        result.obs(i_result)   = scan(isc).obs(iobs).obs; % Computed delay
        result.mjd(i_result)   = mjd; % Delay reference time 
        result.zd_stat1(i_result)    = scan(isc).stat(stat_1_id).zd;
        result.zd_stat2(i_result)    = scan(isc).stat(stat_2_id).zd;
        
%         result.abs_sc_pos_crf(i_result) = norm(sc_pos_crf);
%         result.sc_pos_crf(i_result, :)   = sc_pos_crf';
%         results.vac_del(i_result)        = du;
% 		  result.obs(ierg)   = scan(isc).obs(iobs).obs;
%         results.tpd_g(i_result)        = tpd_g;
%         results.c_trop(i_result)       = c_trop;
%         results.c_therm(i_result)      = c_therm;
%         results.c_axis(i_result)       = c_axis;
    end
        
    end % iobs
    
    scan(isc).space.source =        rq;         % source vector
    scan(isc).space.xp     =        xp;         % pole coordinate x [rad]
    scan(isc).space.yp     =        yp;         % pole coordinate y [rad]
    scan(isc).space.era    =        ERA(isc);   % Earth rotation angle [rad]
    scan(isc).space.xnut   =        XNUT(isc);  % celestial pole X (Xnut+DX)[rad]
    scan(isc).space.ynut   =        YNUT(isc);  % celestial pole Y (Ynut+DY)[rad]
    scan(isc).space.t2c    =        t2c;        % matrix terr2celestial (Q)
    
    % save pantd
    for is=1:length(scan(isc).stat)
        scan(isc).stat(is).pantd = scan(isc).obs(1).pstat2';
    end
    
    % Display counter in CW:
    if mod(isc,100)==0
        disp(['processing scan ', num2str(isc),' of ',num2str(number_of_all_scans)])
    end
    
end % isc


for i_stat=1:length(antenna)
    % save a priori gradients into antenna structure
    antenna(i_stat).apriori_ngr = a_ngr(i_stat)*1000; %[mm]
    antenna(i_stat).apriori_egr = a_egr(i_stat)*1000; %[mm]
    
    
    % #### Error flag messages for eachs station ####
    if flagmess.cto(i_stat) ==1
        fprintf('\n Problems with tidal ocean loading at station %8s \n',antenna(i_stat).name);
    end
    if flagmess.cta(i_stat) ==1
        fprintf('\n Problems with tidal atmosphere loading at station %8s \n',antenna(i_stat).name);
    end
    if flagmess.cnta(i_stat) ==1
        fprintf('\n Problems with non-tidal atmosphere loading at station %8s \n',antenna(i_stat).name);
    end
    
    if flagmess.crg(i_stat) ==1
        fprintf('\n Problems with regression coefficients for atmosphere loading at station %8s \n',antenna(i_stat).name);
    end
    if flagmess.gia(i_stat) ==1
        fprintf('\n Problems with GIA uplift rates at station %8s \n',antenna(i_stat).name);
    end
    
    if flagmess.chl(i_stat) ==1
        fprintf('\n Problems with hydrology loading at station %8s \n',antenna(i_stat).name);
    end
    if flagmess.axtyp(i_stat) ==1
        fprintf('\n No info about mounting type at station %8s. Axis offset correction not applied. \n',antenna(i_stat).name);
    end
    if flagmess.thermal(i_stat) ==1
        fprintf('\n Problems with thermal correction at station %8s \n',antenna(i_stat).name);
    end
    if flagmess.vmf3(i_stat) ==1
        fprintf('\n Problems with VMF3 at station %8s \n',antenna(i_stat).name);
    end
    if flagmess.vmf1(i_stat) ==1
        fprintf('\n Problems with VMF1 at station %8s \n',antenna(i_stat).name);
    end
    if flagmess.dao(i_stat) ==1
        fprintf('\n A priori gradients (DAO) for station %8s are zero\n',antenna(i_stat).name);
    end
    if flagmess.ctop(i_stat) ==1
        fprintf('\n Problems with ocean pole tide loading at station %8s \n',antenna(i_stat).name);
    end
end

% #### Error msg. (station independent) ####
if flgm_ctp == 1
    fprintf('Cubic model after 2010.0 is not available, a linear model for extrapolation is used. (IERS Conv. 2010)\n')
end

% Source structure
if parameter.vie_mod.write_jet==1
    fclose(fidjet);
    fclose(fidjetuv);
    fclose(fidjetjb);
    disp('close fidjet')
end

disp(' ')
disp('vie_mod successfully finished!')


% #################################################
% #### Special handling  and debugging section ####
% #################################################

if strcmp(parameter.vie_mod.special_handling_tag, 'geocentricAprioriDelay')
    
    fprintf(1, '\n');
    fprintf(1, '#############################################################################\n');
    fprintf(1, ' Special handling tag: geocentricAprioriDelay\n');
    fprintf(1, ' => Hit any key to print the geocentric a priori delays for the current session (%s) to the CW!\n',  parameter.session_name);
    fprintf(1, ' => Further instructions:\n');
    fprintf(1, '   1.) Copy all values (the whole column) and append the values as a new column in the input VSO file (%s%s)\n', parameter.filepath, parameter.session_name);
    fprintf(1, '   2.) Use the VSO file to modify the .im files originaly created with calc in the DiFX processing (function vso2im.m in /OUT/).\n');
    fprintf(1, '#############################################################################\n');
    fprintf(1, '  - Number of lines (observations): %d => Make sure that the MATLAB CW is able to show all lines simultaneously\n!', length(result.com));
    fprintf(1, '#############################################################################\n');
    fprintf(1, '\n');
    
    input('Hit any key to continue and write computed geocentric a priori delays to CW!') % Stop here!
    
    format long
    clc
    disp(result.com');
    
    input('Hit any key to continue!') % Stop here!
    
    figure
    plot(result.mjd, result.com*1e9, 'go-')
    title(['computed delay - ', parameter.session_name])
    ylabel('ns')

    % figure
    % plot(result.mjd(1:end-1), diff(result.com*1e9), 'go-')
    % title(['computed delay - ', parameter.session_name])
    % ylabel('ns')
    % 
    % figure
    % plot(result.abs_sc_pos_crf, '-*')
    % 
    % figure
    % plot(-result.com(1:600), '-*')
    % 
    figure
    plot(result.mjd, '-x')
    title('mjd - obs. epochs (modjuldat.m)')
    % 
    % figure
    % plot(sources.s.mjd, '-x')
    % title('mjd - ephem. epochs')
    % 
    % figure
    % plot(result.mjd, result.mjd, '-o')
    % 
    % figure
    % plot(result.mjd-result.mjd(1)*1440, result.com*1e9, 'go-')
    % title(['computed delay - ', parameter.session_name])
    % ylabel('ns')
    %  
    % figure
    % plot(result.mjd-result.mjd(1)*1440, result.obs*1e9, 'bo-')
    % title('observed')
    % ylabel('ns')
    % 
    % figure
    % plot((result.mjd-result.mjd(1))*1440,(result.obs-result.com)*1e9,'ro-')
    % title('observed-computed')
    % xlabel('min')
    % ylabel('ns')
end


