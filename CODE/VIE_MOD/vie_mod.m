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
%       changesa
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
%   25 Nov 2021 by H. Wolf: created some external functions 
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

% Select delay model for quasars:
%   1 = Consensus model (standard in VLBI)
%   2 = Sovers / influence of planets not included yet
%   3 = SSB acceleration (calculation of the time delay following the consensus model + Titov 2011)
delModQ = 1;

% Select delay model for space crafts:
%   1 = Light time equation (geocentric model)
delModS = 1;

% Options for iteration of sat. velocity and position:
ddt_threshold           = 1e-12;    % Iteration ends, if the change in "dt" is less than this value [sec]
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

% ##########################
%       0 Prepare data
% ##########################

constants
global c omega

% Get session name:
session = parameter.session_name;

if parameter.vie_mod.write_jet==1
    fidjet = fopen(parameter.vie_init.jetfilnam,'w');
    fidjetuv = fopen(parameter.vie_init.jetfilnamuv,'w');
    fidjetjb = fopen(parameter.vie_init.jetfilnamjb,'w');
    fprintf('writing jetangles to file\n');
end

% ##### Evaluate the GPT3 to get meteorological data #####
% => Update the scan and antenna structure with met. data (and mapping function parameters from GPT3)
[antenna, scan] = calc_met_data(parameter, antenna, scan);

% ##### CW Output in case met. data from backup models was used #####
% Temp.:
if strcmp(parameter.vie_init.tp,'in situ')  % this if statement is only used in terms of speed, calculation would be the same without
    out_str = {};
    for iStat = 1 : length(antenna)
        if antenna(iStat).gpt3temp
            out_str{length(out_str) + 1} = antenna(iStat).name;
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
    for iStat = 1 : length(antenna)
        if antenna(iStat).gpt3pres
            out_str{length(out_str) + 1} = antenna(iStat).name;
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
    for iStat = 1 : length(antenna)
        if antenna(iStat).gpt3e
            out_str{length(out_str) + 1} = antenna(iStat).name;
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

obsTypeQidx = strcmp({scan.obs_type}, 'q');
obsTypeSidx = strcmp({scan.obs_type}, 's');

if sum(obsTypeQidx) ~= 0
    RA2000(obsTypeQidx)          = deal([sources.q([scan(obsTypeQidx).iso]).ra2000]);
    DE2000(obsTypeQidx)          = deal([sources.q([scan(obsTypeQidx).iso]).de2000]);
    if strncmpi('icrf3',parameter.vie_init.crf(2),5) && (delModQ ~= 3)
        [DE2000, RA2000] = correct_GA(DE2000,RA2000,mean([scan(:).mjd]));
        fprintf(1, 'ICRF3 is used --> GA will be corrected to 2015 using 5.8 muas/year\n');
    end
    sourceNames(obsTypeQidx)     = deal({sources.q([scan(obsTypeQidx).iso]).name});
end
if sum(obsTypeSidx) ~= 0
    sourceNames(obsTypeSidx)     = deal({sources.s([scan(obsTypeSidx).iso]).name});
end

% Prepare time vectors (observation epochs):
MJD         = ([scan.mjd])';                % time of observations MJD UT
LEAP        = tai_utc(MJD);                 % get difference TAI-UTC [s]
TT          = MJD + (32.184 + LEAP)./86400; % [MJD TT]

% ##############################
% ###  1. EARTH ORIENTATION  ###
% ##############################

% -------------------------------------------------
%   READING EARTH ORIENTATION PARAMETERS
% -------------------------------------------------

% Get EOP values for observation epochs:
[MJDeop, XPeop, YPeop, UT1eop, dXeop, dYeop] = load_eop(MJD, parameter); % [d,rad,sec]

% store a priori values at EOP timesteps
parameter.eop.mjd =  MJDeop;
parameter.eop.xp  =  rad2as(XPeop); % [as]
parameter.eop.yp  =  rad2as(YPeop); % [as]
parameter.eop.ut1 =  UT1eop;        % [sec]
parameter.eop.dX  =  rad2as(dXeop); % [as]
parameter.eop.dY  =  rad2as(dYeop); % [as]

% Get EOP values for orbit data epochs:
if ~isempty(sources.s)
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
[parameter,DUT1, XP, YP, DX, DY] = interpolateEOP(parameter, MJDeop, UT1eop, XPeop, YPeop, dXeop, dYeop, MJD, TT, 'observation');
disp('...done.')

% ##### Interpolate EOP values for orbit data epochs: #####

if ~isempty(sources.s) 
    disp('Interpolate EOP values for satellite orbit epochs:')
    [parameter, DUT1_s,XP_s, YP_s,DX_s, DY_s] = interpolateEOP(parameter,MJDeop_s,UT1eop_s,XPeop_s,YPeop_s,dXeop_s,dYeop_s,MJD_s, TT_s, 'orbit');
    disp('...done.')
end


% ##### if ray-tracing files are to be used #####
%set the remaining parameters to VMF3 which can then be used as default, 
%in case a certain line is missing in the .trp-file
raytr_used = strcmp(parameter.vie_init.tropSource.name,'raytr');
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
end 

% --------------------------------------------------------------
%  source vector and derivatives in the catalogue system
% --------------------------------------------------------------

% #### Quasars: ####
% Output: [source vector, partial derivative w.r.t. ra, partial derivative w.r.t. de]
[RQ, DRQDRA, DRQDDE] = sourcevec(DE2000,RA2000); 

% #### Space Crafts: ####
% get crf position and velocity for satellites
if ~isempty(sources.s)
    [sources] = getCRF_PosVelSatellites(sources,T2C_s);
end
% #### Source structure ####
% choose simulated structure catalogue depending on options
if parameter.vie_mod.ssou==1 || parameter.vie_mod.write_jet==1 % have structure
    % read in catalogue
    filename_sou_cat = strcat('../CRF/SOURCE_STRUCTURE_CAT/',parameter.vie_mod.sou_cat);
    [cat_comp.name,cat_comp.flux,cat_comp.maj,cat_comp.min,cat_comp.angle,cat_comp.dRA,cat_comp.dDec] = textread(filename_sou_cat,'%s%f%f%f%f%f%f','delimiter',',');
end

% ##################################
% ###  2. STATION COORDINATES   ####
% ##################################
zz = zeros(length(antenna),1);
flagmess = struct('cto',zz,'cta',zz,'cnta',zz,'crg',zz,'gia',zz,'axtyp',zz,'thermal',zz,'vmf3',zz,'vmf1',zz,'dao',zz,'ctop',zz,'chl',zz);

% antenna corrections for all stations
time_lim = [min(MJD) max(MJD)];     % scans might not be in straight order
[antenna, flagmess] = corr_ant(time_lim, scan(1).tim, antenna, parameter, flagmess);

% antenna positions and velocities
ANT = [[antenna.x]', [antenna.y]', [antenna.z]'];
VEL = [[antenna.vx]', [antenna.vy]', [antenna.vz]'];

[year, month, day, ~, ~, ~] = mjd2date(MJD(1));

% define relevant GRAD data, if needed
if strcmpi(parameter.vie_mod.apgm_h,'grad') || strcmpi(parameter.vie_mod.apgm_w,'grad')
    % read respective file
    gradPath = '../TRP/GRAD/';
    gradFile = [gradPath, 'y', num2str(year), '.grad_r'];
    fidGrad = fopen(gradFile, 'r');
    if fidGrad ~= -1
        gradData = textscan(fidGrad, '%8s  %8.2f  %6.3f %6.3f %6.3f %6.3f', 'CommentStyle', '#','delimiter', '||');
        fclose(fidGrad);
    else
        fprintf('The required GRAD file is not available, no a priori gradients used instead\n');
        parameter.vie_mod.apgm_h = 'no';
        parameter.vie_mod.apgm_w = 'no';
    end
end

% Elliptical coordinates of antennas
[PHI,LAM,H_ELL] = xyz2ell(ANT);

% *************************
%  loop over scans
% *************************
disp('station corrections')

% do that before the loop - otherwise very likely pretty slow!
cpsd_all = cPostSeismDeform(MJD,antenna); % [3 x nScans x nStat] matrix / meters!

if strcmpi(parameter.vie_init.zhd,'gpt3')  || strcmpi(parameter.vie_init.zwd,'gpt3') || strcmpi(parameter.vie_mod.mfh,'gpt3') || strcmpi(parameter.vie_mod.mfw,'gpt3') ||  strcmpi(parameter.vie_mod.apgm_h,'gpt3') || strcmpi(parameter.vie_mod.apgm_w,'gpt3')
    cell_grid_GPT3 = gpt3_5_fast_readGrid;
else 
    cell_grid_GPT3 = 0;
end

% GRAVITATIONAL DEFORMATION
% check if the antenna struct has a field for gravitational deformation. It
% might have been created with an older version of vie_init. If not warn the user.
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

for iSc = 1:number_of_all_scans   
    % running variables for active scan
    mjd   = MJD(iSc);
    tim   = scan(iSc).tim;
    leap  = LEAP(iSc);
    t2c   = T2C(:,:,iSc);
    dQdxp = DQDXP(:,:,iSc);
    dQdyp = DQDYP(:,:,iSc);
    dQdut = DQDUT(:,:,iSc);
    dQddX = DQDDX(:,:,iSc);
    dQddY = DQDDY(:,:,iSc);
    
    sec_of_day = scan(iSc).tim(4)*3600 + scan(iSc).tim(5)*60 + scan(iSc).tim(6);
    
    % tidal ERP variations
    if opt.est_erpsp==1
        PNn     = PNN(:,:,iSc);
        Rall    = RN(:,:,iSc);
        Wn      = WN(:,:,iSc);
        R3ss    = R3SS(:,:,iSc);
        R2xp    = R2XP(:,:,iSc);
        dR2xp   = DR2XP(:,:,iSc);
        R1yp    = R1YP(:,:,iSc);
        dR1yp   = DR1YP(:,:,iSc);
        dRdutn  = DRDUTN(:,:,iSc);
        cosag2  = COSAG2(:,iSc)';
        sinag2  = SINAG2(:,iSc)';
        cosag   = [COSAG1(:,iSc)' COSAG2(:,iSc)'];
        sinag   = [SINAG1(:,iSc)' SINAG2(:,iSc)'];
    end
    
    % FCN partials from nutation
    if opt.est_FCNnut == 1
        dQdfcn = DQDFCN(:,:,iSc);
    end
    if opt.est_FCNnutAmpl == 1
        dQdfcnAc = DQDFCNAC(:,:,iSc);
        dQdfcnAs =DQDFCNAS(:,:,iSc);
    end
    
    % Source:
    % RQ, DRQDRA, DRQDDE only contain valid values for quasar observations.
    % => For spacecraft obs. the values are = 0!
    rq      = RQ(iSc,:);         % Source vector
    drqdra  = DRQDRA(iSc,:);     % Partial derivative w.r.t. Ra
    drqdde  = DRQDDE(iSc,:);     % Partial derivative w.r.t. Dec
    
    % Pole coordinates:
    xp      = XP(iSc);
    yp      = YP(iSc);
    
    % Ephemerides:
    vearth  = ephem.earth(iSc).vbar;
   
    % Compute frequency and phase of 342 tidal constituents, which are
    % recommended for implementation of ocean tidal loading in IERS
    % Conventions 2010.
    [cto_F, cto_P, cto_TAMP, cto_IDD1] = libiers_tdfrph_call(mjd, leap);

    % reference epoch for SSB acceleration
    refep_accSSB = 57023; % 2015.0
    delt_accSSB  = (mjd - refep_accSSB) * 86400; % time since ref epoch in sec
 
    % ************************************
    %  loop over stations in current scan
    % ************************************
    for iStat = 1 : length(scan(iSc).stat)
       [scan, flgm_ctp] = antennaCorrections(iSc, iStat, scan, antenna, opt, parameter, session, mjd, ANT, VEL, t2c, leap, cto_F, cto_P, cto_TAMP, cto_IDD1, PHI, LAM, tim, xp, yp, cpsd_all, ephem);    
    end 
    
    % ############################
    % ###  3. COMPUTED DELAY   ###
    % ############################
       
    % ****************************************************
    %  loop over baselines (observations) in current scan
    % ****************************************************
    for iobs = 1:length(scan(iSc).obs)
        % Partial derivatives of the delay w.r.t. the position of space crafts (different ref. frames)           
        pdSatPosGCRF  = [];
        pdSatPosTRF   = [];
        pdSatPosRSW   = [];
        
        % get station IDs of the current baseline:
        idStation1  = scan(iSc).obs(iobs).i1;
        idStation2  = scan(iSc).obs(iobs).i2;
        
        % get station information
        crsStation1 = scan(iSc).stat(idStation1).xcrs;    % station1 position CRS
        crsStation2 = scan(iSc).stat(idStation2).xcrs;    % station2 position CRS
        trsStation1  = scan(iSc).stat(idStation1).x;       % station1 position TRS
        trsStation2  = scan(iSc).stat(idStation2).x;       % station2 position TRS
        
        rqu = rq / norm(rq);                            % unit source vector barycentrum-source, only valid for quasar obs.! For spacecrafts = [1,0,0]
        crsBaseline  = crsStation2 - crsStation1;       % baseline vector CRS
        trsBaseline  = trsStation2-trsStation1;         % baseline vector TRS
        
        % station velocity due to earth rotation
        v1 = t2c *[-omega*trsStation1(2);omega*trsStation1(1);0];  % [CRS]
        v2 = t2c *[-omega*trsStation2(2);omega*trsStation2(1);0];  % [CRS]
             
        % ##### Distinguish between source types (quasar/spacecraft) #####
        delModQ = 1;
        switch(scan(iSc).obs_type)
            case 'q'         
                switch delModQ
                    case 1 %consensus        
                        [tau, pGammaSun, k1a, k2a, fac1] = consensusModelQuasar(iSc, crsStation1,crsStation2, idStation1, idStation2, rqu, ephem, opt, antenna, v1, v2);

                    case 2 % sovers / influence of planets not included yet
                        [tau, pGammaSun, k1a, k2a, fac1] = soversModelQuasar(iSc, crsStation1,crsStation2, rqu, ephem, v2);
                        
                    case 3 % SSB acceleration  
                        [tau, pGammaSun, k1a, k2a, fac1] = SSBAccelerationModelQuasar(iSc, crsStation1,crsStation2, rqu, ephem, v1, v2, opt, delt_accSSB);
                end
            case 's'
                switch(delModS)                
                    case 1 % Light time equation
                        [ps1, ps2, pdSatPosRSW, pdSatPosGCRF, pdSatPosTRF, tau] = calcDelaySatellite(sources, iSc, scan, idStation1, idStation2, flag_fix_sat_pos_to_stat1, ddt_threshold, crsStation1, crsStation2, antenna, sec_of_day, mjd, max_iterations, ephem, v2, t2c);
                end        
        end 
        
        % further corrections (same for both models (Sekido & Fukushima, p.141))
        [a_ngr, a_egr, scan, antenna, tau]  = correctionBaseline(scan, antenna, parameter, t2c, mjd, iSc, idStation1, idStation2, k1a, k2a, rqu, v2, v1, tau, cell_grid_GPT3, session, iobs);
            
        % SOURCE STRUCTURE +
        if parameter.vie_mod.ssou==1 || parameter.vie_mod.write_jet==1
            ind=strcmp(sources.q(scan(iSc).iso).name,cat_comp.name,'exact');
            if isempty(ind)
                nam2=sounamivs2iers(sources.q(scan(iSc).iso).name);
                ind = strcmp(nam2,cat_comp.name,'exact');
            end
            if isempty(ind)
                disp(strcat('source ',sources.q(scan(iSc).iso).name,' not found in ss catalogue.'));
                soucorr=0; jetang=90; jetjb=0; uvrange=0; uu=0; vv=0;
            else
                sources.q(scan(iSc).iso).sou_model=[cat_comp.flux(ind),cat_comp.maj(ind),cat_comp.min(ind),cat_comp.angle(ind),cat_comp.dRA(ind),cat_comp.dDec(ind)];
                [soucorr,uu,vv]=modDelay(sources.q(scan(iSc).iso).sou_model,scan(iSc).stat(idStation1).x,scan(iSc).stat(idStation2).x,([8213 8252 8353 8513 8733 8853 8913 8933]+4), 8217, sources.q(scan(iSc).iso).ra2000, sources.q(scan(iSc).iso).de2000, scan(iSc).mjd);
                soucorr=soucorr*1e-12;
                jetvec=(sources.q(scan(iSc).iso).sou_model(2,5:6));
                jetvec=jetvec/norm(jetvec);
                uvvec=[uu;vv];
                uvvec=uvvec/norm(uvvec);
                jetang=acos(abs(jetvec*uvvec))*180/pi;
                uvrange=(jetvec)*[uu;vv];
                jetjb=(90-jetang)*(norm(crsBaseline)/6371000);
            end
        end
        
        % correct for source structure
        if parameter.vie_mod.ssou==1
            tau=tau+soucorr;
        end
        
        % writing jetang to external file
        if parameter.vie_mod.write_jet==1
            fprintf(fidjet,'%s %s %4.12f %4.12f\n',antenna(idStation1).name,antenna(idStation2).name,mjd,jetang);
            fprintf(fidjetuv,'%s %s %4.12f %4.12f\n',antenna(idStation1).name,antenna(idStation2).name,mjd,uvrange);
            fprintf(fidjetjb,'%s %s %4.12f %4.12f\n',antenna(idStation1).name,antenna(idStation2).name,mjd,jetjb);
        end
        % SOURCE STRUCTURE -         
        
        % ##############################
        % ### 4. PARTIAL DERIVATIVES ###
        % ##############################
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
        M   = (dij -(rq'*b2')/rho)*(-gam*(1-b2'*beta)*(E*crsBaseline)/rho); %[3,1]
        
        % wrt EOP per baseline [m]
        pdx = K'* (dQdxp * trsBaseline')/c;
        pdy = K'* (dQdyp * trsBaseline')/c;
        put = K'* (dQdut * trsBaseline')/c;
        pdX = K'* (dQddX * trsBaseline')/c;
        pdY = K'* (dQddY * trsBaseline')/c;
        
        % wrt source coordinates [cm/mas]
        % RaDec position of quasars
        % - set=0, if no quasar was observed in this scan (e.g. scan to spacecraft)
        if strcmp(scan(iSc).obs_type, 'q')
            parra = (drqdra * M)*pi()/180/3600000*100;  % (2.230)
            parde = (drqdde * M)*pi()/180/3600000*100;  % (2.231)  % [cm/mas]
			psou = [parra, parde];
        else 
            psou = [];
        end
        
        % wrt station coordinates
        if strcmp(scan(iSc).obs_type, 'q')
            ps1  = -B; % (2.248)
            ps2  =  B; % (2.248)
        end
        
        % Axis Offset
        paxktStation1 = -scan(iSc).stat(idStation1).daxkt; %[-]
        paxktStation2 =  scan(iSc).stat(idStation2).daxkt; %[-]
 
        % Amplitudes of the seasonal variation of station position
        pAcrStation1 = (- scan(iSc).stat(idStation1).pAcr_xyz)* B';
        pAceStation1 = (- scan(iSc).stat(idStation1).pAce_xyz)* B';
        pAcnStation1 = (- scan(iSc).stat(idStation1).pAcn_xyz)* B';
        pAsrStation1 = (- scan(iSc).stat(idStation1).pAsr_xyz)* B';
        pAseStation1 = (- scan(iSc).stat(idStation1).pAse_xyz)* B';
        pAsnStation1 = (- scan(iSc).stat(idStation1).pAsn_xyz)* B';
        
        pAcrStation2 = (scan(iSc).stat(idStation2).pAcr_xyz)* B';
        pAceStation2 = (scan(iSc).stat(idStation2).pAce_xyz)* B';
        pAcnStation2 = (scan(iSc).stat(idStation2).pAcn_xyz)* B';
        pAsrStation2 = (scan(iSc).stat(idStation2).pAsr_xyz)* B';
        pAseStation2 = (scan(iSc).stat(idStation2).pAse_xyz)* B';
        pAsnStation2 = (scan(iSc).stat(idStation2).pAsn_xyz)* B';
        
        % pole tide Love & Shida per baseline
        phpole_bl = (scan(iSc).stat(idStation2).phpole - scan(iSc).stat(idStation1).phpole)* B';
        plpole_bl = (scan(iSc).stat(idStation2).plpole - scan(iSc).stat(idStation1).plpole)* B';
          
        % APL RG
        prgStation1 = (-scan(iSc).stat(idStation1).prg)* B';
        prgStation2 = (scan(iSc).stat(idStation2).prg)* B';
        
        % tidal ERP partials [m]     
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
                    R2xp*dR1yp.*sinag(tidn)))) * trsBaseline')/c;
                pbp(tidn) = K'* ((PNn*Rall* (R3ss*(dR2xp.*sinag(tidn)*R1yp + ...
                    R2xp*dR1yp.*cosag(tidn)))) * trsBaseline')/c;
                puc(tidn) = K'* ((PNn    * dRdutn.*cosag(tidn) * Wn) * trsBaseline')/c;
                pus(tidn) = K'* ((PNn    * dRdutn.*sinag(tidn) * Wn) * trsBaseline')/c;
            end
            
            % retrograde semidiurnal polar motion
            for tidn = 1 : length(cosag2)
                pam(tidn) = K'* ((PNn*Rall* (R3ss*(dR2xp.*(-cosag2(tidn))*R1yp + ...
                    R2xp*dR1yp.*(-sinag2(tidn))))) * trsBaseline')/c;
                pbm(tidn) = K'* ((PNn*Rall* (R3ss*(dR2xp.*(-sinag2(tidn))*R1yp + ...
                    R2xp*dR1yp.*cosag2(tidn)))) * trsBaseline')/c;
            end
        end
        
        if opt.est_FCNnut == 1
            pFCNnut = K'* (dQdfcn * trsBaseline')/c;
        end
        if opt.est_FCNnutAmpl == 1
            pFCNnutAc = K'* (dQdfcnAc * trsBaseline')/c; %[s]
            pFCNnutAs = K'* (dQdfcnAs * trsBaseline')/c; %[s]
        end
        
        % Love & Shida per baseline
        pLove_bl  = (scan(iSc).stat(idStation2).pLove  - scan(iSc).stat(idStation1).pLove)* B';
        pShida_bl = (scan(iSc).stat(idStation2).pShida - scan(iSc).stat(idStation1).pShida)* B';
        pFCN_bl   = (scan(iSc).stat(idStation2).pFCN   - scan(iSc).stat(idStation1).pFCN)* B';
        
        % wrt SBB acceleration
        pacc = delt_accSSB/c^2*((rq*crsBaseline)*rq'-crsBaseline)...
            -delt_accSSB/c^3*((rq*(vearth+v2))*crsBaseline...
            +((rq*vearth)*crsBaseline)/2 - (crsBaseline'*vearth)*rq'); %[sec^3/m]
        pacc = pacc'./100; %[sec^3/cm]
        
        % ############################
        % ###  5. STORE RESULTS   #### 
        % ############################
        
        scan(iSc).obs(iobs).com = tau; %[sec]
        
        % partial derivatives of the delay  w.r.t.: 
        scan(iSc).obs(iobs).psou                = psou;             % source coordinates [cm/mas]; =[0,0], if no quasar was observed
        scan(iSc).obs(iobs).psou_sat_gcrf       = pdSatPosGCRF;     % satellite coordinates in GCRF [sec/m]
        scan(iSc).obs(iobs).psou_sat_trf        = pdSatPosTRF;      % satellite coordinates in TRF [sec/m]
        scan(iSc).obs(iobs).psou_sat_rsw        = pdSatPosRSW;      % satellite coordinates in RSW system ("satellite coord. sys.") [sec/m]
        scan(iSc).obs(iobs).pnut                = [pdX,pdY];        % dX, dY [sec/rad]
        scan(iSc).obs(iobs).ppol                = [pdx,pdy,put];    % xpol, ypol, dut1 [sec/rad]
        scan(iSc).obs(iobs).pstat1              = ps1;              % station1
        scan(iSc).obs(iobs).pstat2              = ps2;              % station2
        scan(iSc).obs(iobs).pAO_st1             = paxktStation1;    % Axis offset at station 1 [-]
        scan(iSc).obs(iobs).pAO_st2             = paxktStation2;    % Axis offset at station 2 [-]
       
        scan(iSc).obs(iobs).pAcr_st1 = pAcrStation1;  % Amplitudes of the seasonal variation in station pos.
        scan(iSc).obs(iobs).pAce_st1 = pAceStation1;  % [-]
        scan(iSc).obs(iobs).pAcn_st1 = pAcnStation1;
        scan(iSc).obs(iobs).pAsr_st1 = pAsrStation1;
        scan(iSc).obs(iobs).pAse_st1 = pAseStation1;
        scan(iSc).obs(iobs).pAsn_st1 = pAsnStation1;
        
        scan(iSc).obs(iobs).pAcr_st2 = pAcrStation2;  % Amplitudes of the seasonal variation in station pos.
        scan(iSc).obs(iobs).pAce_st2 = pAceStation2;  % [-]
        scan(iSc).obs(iobs).pAcn_st2 = pAcnStation2;
        scan(iSc).obs(iobs).pAsr_st2 = pAsrStation2;
        scan(iSc).obs(iobs).pAse_st2 = pAseStation2;
        scan(iSc).obs(iobs).pAsn_st2 = pAsnStation2;
        
        scan(iSc).obs(iobs).phpole   = phpole_bl; % [cm] for pole tide Love number
        scan(iSc).obs(iobs).plpole   = plpole_bl; % [cm] for pole tide Shida number
        
        scan(iSc).obs(iobs).prg_st1  = prgStation1; %[hPa]  APL RG
        scan(iSc).obs(iobs).prg_st2  = prgStation2; %[hPa]	APL RG
        
        % tidal ERP
        if opt.est_erpsp == 1
            scan(iSc).obs(iobs).ppap  = pap;  % prograde pm cos
            scan(iSc).obs(iobs).ppam  = pam;  % retrograde pm cos
            scan(iSc).obs(iobs).ppbp  = pbp;  % prograde pm sin
            scan(iSc).obs(iobs).ppbm  = pbm;  % retrograde pm sin
            scan(iSc).obs(iobs).putc  = puc;  % ut1 cos
            scan(iSc).obs(iobs).puts  = pus;  % ut1 sin
        end
        
        if opt.est_FCNnut ==1
            scan(iSc).obs(iobs).pFCNnut     =  pFCNnut;     % [sec*day]
        else
            scan(iSc).obs(iobs).pFCNnut     =  0;           % [sec*day]
        end
        if opt.est_FCNnutAmpl==1
            scan(iSc).obs(iobs).pFCNnutAc   =  pFCNnutAc;   % [sec]
            scan(iSc).obs(iobs).pFCNnutAs   =  pFCNnutAs;   % [sec]
        else
            scan(iSc).obs(iobs).pFCNnutAc   =  0;           % [sec]
            scan(iSc).obs(iobs).pFCNnutAs   =  0;           % [sec]
        end
        
        scan(iSc).obs(iobs).pLove    = pLove_bl;    % [cm] for Love numbers
        scan(iSc).obs(iobs).pShida   = pShida_bl;   % [cm] for Shida
        scan(iSc).obs(iobs).pFCN     = pFCN_bl;     % [cm] for FCN
        scan(iSc).obs(iobs).pacc     = pacc;        % dt/dacc SSB acceleration [sec^3/cm]
        scan(iSc).obs(iobs).pGamma   = pGammaSun;   % [sec]       
        scan(iSc).obs(iobs).pscale   = -fac1;       % [sec] Correction to the scale factor

        if flag_save_results
            i_result = i_result + 1;
            result.com(i_result)   = tau; % Computed delay
            result.obs(i_result)   = scan(iSc).obs(iobs).obs; % Computed delay
            result.mjd(i_result)   = mjd; % Delay reference time 
            result.zd_stat1(i_result)    = scan(iSc).stat(idStation1).zd;
            result.zd_stat2(i_result)    = scan(iSc).stat(idStation2).zd;

    %         result.abs_sc_pos_crf(i_result) = norm(sc_pos_crf);
    %         result.sc_pos_crf(i_result, :)   = sc_pos_crf';
    %         results.vac_del(i_result)        = du;
    % 		  result.obs(ierg)   = scan(iSc).obs(iobs).obs;
    %         results.tpd_g(i_result)        = tpd_g;
    %         results.c_trop(i_result)       = c_trop;
    %         results.c_therm(i_result)      = c_therm;
    %         results.c_axis(i_result)       = c_axis;
        end
        
    end 
    scan(iSc).space.source =        rq;         % source vector
    scan(iSc).space.xp     =        xp;         % pole coordinate x [rad]
    scan(iSc).space.yp     =        yp;         % pole coordinate y [rad]
    scan(iSc).space.era    =        ERA(iSc);   % Earth rotation angle [rad]
    scan(iSc).space.xnut   =        XNUT(iSc);  % celestial pole X (Xnut+DX)[rad]
    scan(iSc).space.ynut   =        YNUT(iSc);  % celestial pole Y (Ynut+DY)[rad]
    scan(iSc).space.t2c    =        t2c;        % matrix terr2celestial (Q)
    
    % save pantd
    for is=1:length(scan(iSc).stat)
        scan(iSc).stat(is).pantd = scan(iSc).obs(1).pstat2';
    end
    
    % Display counter in CW:
    if mod(iSc,100)==0
        disp(['processing scan ', num2str(iSc),' of ',num2str(number_of_all_scans)])
    end
end

for iStat=1:length(antenna)
    % save a priori gradients into antenna structure
    antenna(iStat).apriori_ngr = a_ngr(iStat)*1000; %[mm]
    antenna(iStat).apriori_egr = a_egr(iStat)*1000; %[mm]
    
    % #### Error flag messages for each station ####
    if flagmess.cto(iStat) ==1
        fprintf('\n Problems with tidal ocean loading at station %8s \n',antenna(iStat).name);
    end
    if flagmess.cta(iStat) ==1
        fprintf('\n Problems with tidal atmosphere loading at station %8s \n',antenna(iStat).name);
    end
    if flagmess.cnta(iStat) ==1
        fprintf('\n Problems with non-tidal atmosphere loading at station %8s \n',antenna(iStat).name);
    end
    
    if flagmess.crg(iStat) ==1
        fprintf('\n Problems with regression coefficients for atmosphere loading at station %8s \n',antenna(iStat).name);
    end
    if flagmess.gia(iStat) ==1
        fprintf('\n Problems with GIA uplift rates at station %8s \n',antenna(iStat).name);
    end
    
    if flagmess.chl(iStat) ==1
        fprintf('\n Problems with hydrology loading at station %8s \n',antenna(iStat).name);
    end
    if flagmess.axtyp(iStat) ==1
        fprintf('\n No info about mounting type at station %8s. Axis offset correction not applied. \n',antenna(iStat).name);
    end
    if flagmess.thermal(iStat) ==1
        fprintf('\n Problems with thermal correction at station %8s \n',antenna(iStat).name);
    end
    if flagmess.vmf3(iStat) ==1
        fprintf('\n Problems with VMF3 at station %8s \n',antenna(iStat).name);
    end
    if flagmess.vmf1(iStat) ==1
        fprintf('\n Problems with VMF1 at station %8s \n',antenna(iStat).name);
    end
    if flagmess.dao(iStat) ==1
        fprintf('\n A priori gradients (DAO) for station %8s are zero\n',antenna(iStat).name);
    end
    if flagmess.ctop(iStat) ==1
        fprintf('\n Problems with ocean pole tide loading at station %8s \n',antenna(iStat).name);
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


