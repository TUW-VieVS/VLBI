% ************************************************************************
%   Description:
%   function to create a default parameterization - structure array
%   'opt' is read from 'multi_lsmopt.m' in 'VIE_SETUP_V1a'
%
%   Reference: 
%
%   Input:	
%       'antenna'           structure array     (for info. /DOC/antenna.doc)
%       'sources'           structure array     (for info. /DOC/sources.doc)
%       'na'                (1,1)               number of antennas
%       'ns'                (1,1)               number of sources
%       'obs_per_source'    structure array     (for info. /DOC/obs_per_source.doc)
%       'opt'        structure array     (for info. /DOC/opt.doc)
%
%   Output:
%       'opt'        structure array     (for info. /DOC/opt.doc)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   06 July 2009 by Johannes Boehm
%
%   Revision: 
%   14 Aug 2009 by Kamil Teke: options for global estmimates added
%   12 Sep 2009 by Kamil Teke: options for sources added
%   23 Nov 2009 by Kamil Teke: GUI for options added
%   06 Dec 2009 by Kamil Teke: header added
%   05 Feb 2010 by Tobias Nilsson: Removed small bug occuring when line
%      with reference clock in OPT-file is shorter than 8 characters
%   05 Mar 2010 by Kamil Teke: options of troposphere gradients absolute constraints added
%   12 Oct 2010 by Hana Spicakova: opt.sinex_out for N and b; prepares N
%      and b for future sinex output
%   11 Nov 2010 by Hana Spicakova: FCN period and reference epoch for station
%                             velocities added
%   12 Nov 2010 by Hana Spicakova: indicator for single-session solution added
%                             and opt.acsii_snx for ASCII SINEX output
%   17 May 2011 by Hana Spicakova: opt.sinex_out deleted from the
%                             structure, it is not needed any more
%   17 May 2011 by Hana Spicakova: variables opt.snxout.clk, opt.snxout.zwd,...
%                             added for the SINEX output
%   24 March 2013 by Hana Krasna: IERS name added to opt.source
%   10 Jan 2014 by Hana Krasna: AO and APL RgC added
%   19 Jun 2014 by Hana Krasna: sources can be estimated in vie_lsm with
%   NNR condition
%   25 Jun 2014 by Hana Krasna: bug corrected - in a 2 station session where
%   both stations have a clock break the reference clock is not changed 
%   21 Oct 2015 by Hana Krasna: seasonal variations in station positions (annual, semi-annual), pole
%               tide Love and Shida numbers, APL regression coefficient added
%   01 Mar 2017 by A. Hellerschmied: "opt" is not gloabl any more => added as in/output argument
% ************************************************************************


function opt = lsmopt_scan_global(antenna,na,obs_per_source,ns,sources, opt)
% CLAUDIA 10 sept. 2012
% ----------------------
% --- first solution ---
% ----------------------

% first simple parameterization to remove large clock offsets, rates, and
% squared rates; known clock offsets should be taken into account
% this solution can also be used to detect clock breaks
opt.first = opt.first;  % 1 if a first basic solution will be carried out 
                        % 0 if not 

% if optfirst == 1 then
% which clock parameters should be estimated in the first solution
% one reference clock is deleted here (no NNT for the clocks here)
% a station with a clock break must not be used as reference clock
opt.firstclock = opt.firstclock;     % 0 only one offset per station   
                                     % 1 one offset and one rate per station 
                                     % 2 one offset, one rate and one squared rate per station 
% no other estimates possible in this first step
% the clock values found here are used to correct the observations
% clock breaks can be found if residuals per station are plotted

opt.treat_breaks = opt.treat_breaks; % 0 not take clock breaks in to account in the first solution 
                                     % 1 take clock breaks in to account in the firstsolution 

% ------------------------------
% --- second (main) solution ---
% ------------------------------

opt.second = opt.second; % 1 if main solution is carried out 
                         % 0 if no main solution is carried out 

% ---                        
% clock parameterization
opt.pw_clk = opt.pw_clk;     % 0 no clocks estimated 
                             % 1 only pw linear clock offsets 
                             % 2 one rate per session in addition to pwl clock offsets 
                             % 3 one squared rate per station in addition to (2) 
                        
    % case 1,2,3 constraints between clock offsets
    opt.constr_clk = opt.constr_clk ;     % 1 with constraints 
                                          % 0 no constraints 
                             
% ---                 
% zenith delays
opt.pw_zwd = opt.pw_zwd;         % 0 no zwd estimated 
                                 % 1 zenith wet delays estimated as pwl offsets
    % case 1 constraints between zenith delays
    opt.constr_zwd = opt.constr_zwd;     % 1 with constraints between zenith delays 
                                         % 0 no constraints between zd 
                            
% ---                            
% north gradients
opt.pw_ngr = opt.pw_ngr;         % 0 no ngr estimated 
                                 % 1 ngr estimated as pwl offsets 
    % case 1 constraints between north gradients
opt.constr_rel_ngr = opt.constr_rel_ngr;   % 1 with relative constraints between north gradients 
                                           % 0 no relative constraints between north gradients 
                                         
    
opt.constr_abs_ngr = opt.constr_abs_ngr;  % absolute constraints between north gradients 
                            
% ---                            
% east gradients
opt.pw_egr = opt.pw_egr;         % 0 no egr estimated 
                                 % 1 egr estimated as pwl offsets
    % case 1 constraints between east gradients
opt.constr_rel_egr = opt.constr_rel_egr;    % 1 with relative constraints between east gradients 
                                            % 0 no relative constraints between east gradients 
                                         
    
opt.constr_abs_egr = opt.constr_abs_egr;  % absolute constraints between east gradients 
                            
% ---        
% station coordinates                    
% how to estimate station coordinates 
opt.stc = opt.stc;         % 0 not estimate station coordinates (stationwise) 
                           % 1 estimate station coordinates (stationwise)
opt.pw_stc = opt.pw_stc;      % 0 estimate all selected stations as one offset per session (NNT/NNR is available)
                              % 1 estimate all selected station coordinates as pwl offsets (NNT/NNR is NOT available - fix at least 3 stations)
    % case 1 introduce constraints between station coordinates (only if pwl offsets)
    opt.constr_xyz = opt.constr_xyz;     % 1 constraints between pwl coordinate offsets 
                                         % 0 no constraints between pwl coordinate offsets 
    % case 2 NNT for station coordinates
    opt.nnt_stc = opt.nnt_stc ;        % 1 NNT for stations 
                                       % 0 no NNT for stations 
    % case 2 NNR for station coordinates
    opt.nnr_stc = opt.nnr_stc;        % 1 NNR for stations 
                                      % 0 no NNR for stations 
    % case 2 Scale for station coordinates
    opt.sca_stc = opt.sca_stc;        % 1 Scale for stations 
                                      % 0 no Scale for stations 
                                      
% ---                           
% source coordinates                    
% how to estimate source coordinates 
opt.pw_sou = opt.pw_sou;      % 0 not estimate coordinates of source's (piecewise - single session solution)
                              % 1 estimate source's coordinates ... (piecewise - single session solution) 
                     
    opt.constr_sou = opt.constr_sou;     % 1 constraints between pwl source coordinates offsets 
                                         % 0 no constraints between pwl
                                         % offsets of sources coordinates

% ---                            
% in future there will be the possibility to read the following
% station-dependent information from a station master file
for istat = 1:na
    opt.stat(istat).name = antenna(istat).name;  
                                            % station name   
    opt.stat(istat).zwd_inc = opt.pw_zwd;   % 1 include the zwd parameters of the station 
                                            % 0 do not include 
    opt.stat(istat).ngr_inc = opt.pw_ngr;   % 1 include the north gradients parameters of the station 
                                            % 0 do not include 
    opt.stat(istat).egr_inc = opt.pw_egr;   % 1 include the east gradients parameters of the station 
                                            % 0 do not include 
    opt.stat(istat).xyz_inc = opt.stc;      % 0 not estimate station coordinates (stationwise) 
                                            % 1 estimate station coordinates (stationwise)
    % constraint values if applicable
    opt.stat(istat).coef_clk = opt.coef_clk;                 % clock constraints in ps**2/s 
    opt.stat(istat).coef_zwd = opt.coef_zwd;                 % zwd constraints in ps**2/s 
    opt.stat(istat).coef_rel_ngr = opt.coef_rel_ngr;         % ngr relative constraints in mm/day 
    opt.stat(istat).coef_rel_egr = opt.coef_rel_egr;         % egr realtive constraints in mm/day 
    opt.stat(istat).coef_abs_ngr = opt.coef_abs_ngr;         % ngr absolute constraints in mm 
    opt.stat(istat).coef_abs_egr = opt.coef_abs_egr;         % egr absolute constraints in mm 
    opt.stat(istat).coef_xyz = opt.coef_xyz;                 % coordinate constraints in mm/day 
    % time intervals if applicable
    opt.stat(istat).int_clk = opt.int_clk;          % interval of clock offsets in minutes 
    opt.stat(istat).int_zwd = opt.int_zwd;          % interval of wet zenith delay offsets in minutes 
    opt.stat(istat).int_ngr = opt.int_ngr;          % interval of north gradients offsets in minutes 
    opt.stat(istat).int_egr = opt.int_egr;          % interval of east gradients offsets in minutes 
    opt.stat(istat).int_xyz = opt.int_xyz;          % interval of coordinates in offsets minute 
    % other settings
    opt.stat(istat).ref = opt.ref;          % 1 if reference clock 
                                            % 0 if not reference clock 
    opt.stat(istat).nnt_inc = opt.nnt_stc;  % 1 include station in translation minimum 
                                            % 0 do not include 
    opt.stat(istat).nnr_inc = opt.nnr_stc;  % 1 include station in rotation minimum 
                                            % 0 do not include 
    opt.stat(istat).nns_inc = opt.sca_stc;  % 1 include station in scale minimum 
                                            % 0 do not include 
end

% Boehm 19 August 2009 (Excluding station from NNT/NNR)
for istat = 1:na
    if antenna(istat).in_trf == 0
        opt.stat(istat).nnt_inc = 0;
        opt.stat(istat).nnr_inc = 0;
        opt.stat(istat).nns_inc = 0;
        % station is estimated
	opt.stc = 1;
	opt.pw_stc = 0;
        opt.stat(istat).xyz_inc = 1;  % could also be 1
        disp([antenna(istat).name,' is excluded from datum or is not fixed']);
    end
end

% ---
% global parameters: these are only added for global solution but not for
% session-wise solution
% Global solution 
opt.global_solve = opt.global_solve ; % 1 form and save N matrix, b vector and auxiliary info. 
                                      % for global solution 

% Sinex output
opt.ascii_snx = opt.ascii_snx;     % 1 creates ASCII SINEX output (hana 12Nov2010)

% Parameters in the SINEX file
opt.outsnx.clk = opt.outsnx.clk;   % 1 include in SINEX; 0 reduce from SINEX (hana 17May2011)
opt.outsnx.zwd = opt.outsnx.zwd;
opt.outsnx.tgr = opt.outsnx.tgr;
opt.outsnx.sou = opt.outsnx.sou;
opt.outsnx.xyz = opt.outsnx.xyz;
opt.outsnx.eop = opt.outsnx.eop;

% seasonal variations in station position
opt.est_stsespos=opt.est_stsespos;
% Love number for pole tide
opt.est_hpole=opt.est_hpole;
% Shida number for pole tide
opt.est_hpole=opt.est_hpole;

% APL regression coefficients
opt.est_rg = opt.est_rg;


% Love numbers, set up either all or none
opt.est_love = opt.est_love;       % 1 all Love numbers 
                                   % 0 no Love numbers 
% Shida numbers
opt.est_shida = opt.est_shida;      % 1 all Shida numbers 
                                    % 0 no Shida numbers 
% Free Core Nutation period from Solid Earth tides
opt.est_FCNset = opt.est_FCNset;        % 1/0

% acceleration of SSB
opt.est_accSSB=opt.est_accSSB; 
% velocity of sources
opt.est_source_velo=opt.est_source_velo;


% Axis offset
opt.est_AO=opt.est_AO;


% Source coordinates
opt.est_source = opt.est_source;     % 1 all sources coordinates  
                                     % 0 no source coordinates            
% Station velocities
opt.est_vel = opt.est_vel;  % 1 all station velocities
                            % 0 no station velocities
% Reference epoch for station velocities (hana 11Nov2010)
opt.refvel = opt.refvel;      % reference epoch in years

                            
                            
% ---                            
% in future there will be the possibility to read the following
% source-dependent information from a source master file    
                        
% Kamil 27 August 2009
for isou = 1 :ns
    opt.source(isou).total_obs = length(obs_per_source(isou).iso);
end
opt.limit_num_obs_sou = max([opt.source.total_obs]); % total number of obs. per source
for isou = 1 : ns
    opt.source(isou).name = sources(isou).name; % name of the source
    opt.source(isou).IERSname = sources(isou).IERSname; % IERS name of the source (hana March24 2013)
    opt.source(isou).total_obs = length(obs_per_source(isou).iso);
    opt.source(isou).coef_rade = opt.sour_coef_rade; % source coordinate constraints in mas 
    opt.source(isou).int_rade = opt.sour_int_rade; % estimation intervals [minutes] 
    opt.source(isou).rade_inc = opt.pw_sou;
    opt.source(isou).de = sources(isou).de2000;
    opt.source(isou).ra = sources(isou).ra2000;
    opt.source(isou).nnr_inc = opt.est_sourceNNR; % NNR condition on sources
    if opt.global_solve == 0 | opt.global_solve == 1 | opt.est_sourceNNR == 1 & opt.est_source == 0
        if opt.pw_sou == 1
             %if length([obs_per_source(isou).i1]) == opt.limit_num_obs_sou
             %% estimating the source that have maximum number of
             %observations
             if (sources(isou).in_crf==0)
                 opt.source(isou).rade_inc = opt.pw_sou;
             else 
                 opt.source(isou).rade_inc = 0;
             end
        end
        if opt.global_solve == 1 & opt.est_source == 1
                opt.source(isou).rade_inc = opt.est_source; 
        end
        if opt.est_sourceNNR == 1
            if (sources(isou).in_crf==0)
                opt.source(isou).nnr_inc = 0;
            end
        end
    end
end
    
% ---
% Earth orientation parameters
% XPOL
opt.xpol.model = opt.xpol.model;         % 1 piecewise linear offsets  
                                         % 0 not estimated 
    % case 1                        
    opt.xpol.int = opt.xpol.int;        % case 1 interval in minutes  
    % case 1 
    opt.xpol.constrain = opt.xpol.constrain;     % case 1 if constraints should be applied
                                                 % 0 if not applied
        % case 1 
        opt.xpol.coef = opt.xpol.coef;          % constraint in mas/day 
                                                % (e.g. 0.0001 creates just one offset, 3 is loose)
% YPOL
opt.ypol.model = opt.ypol.model;         % 1 piecewise linear offsets 
                                         % 0 not estimated 
    % case 1                         
    opt.ypol.int = opt.ypol.int;        % case 1 interval in minutes  
    % case 1 
    opt.ypol.constrain = opt.ypol.constrain;     % case 1 if constraints should be applied 
                                                 % 0 if not applied 
        % case 1 
        opt.ypol.coef = opt.ypol.coef;          % constraint in mas/day  
                                                % (e.g. 0.0001 creates just one offset, 3 is loose)
% dUT1
opt.dut1.model = opt.dut1.model;         % 1 piecewise linear offsets  
                                         % 0 not estimated 
    % case 1                        
    opt.dut1.int = opt.dut1.int;        % case 1 interval in minutes  
    % case 1 
    opt.dut1.constrain = opt.dut1.constrain;     % case 1 if constraints should be applied
                                                 % 0 if not applied
        % case 1 
        opt.dut1.coef = opt.dut1.coef;          % constraint in mas/day 
                                                % (e.g. 0.0001 creates just one offset, 3/15 is loose)
% Nutation DEPS (Celestial pole offset, dx) 
opt.nutdx.model = opt.nutdx.model;          % 1 piecewise linear offsets 
                                            % 0 not estimated 
    % case 1                          
    opt.nutdx.int = opt.nutdx.int;         % case 1 interval in minutes  
    % case 1 
    opt.nutdx.constrain = opt.nutdx.constrain;      % case 1 if constraints should be applied 
                                                    % 0 if not applied 
        % case 1 
        opt.nutdx.coef = opt.nutdx.coef;      % constraint in mas/day 
                                     % (e.g. 0.0001 creates just one offset, 3 is loose)
% Nutation DPSI (Celestial pole offset, dy)
opt.nutdy.model = opt.nutdy.model;        % 1 piecewise linear offsets 
                                          % 0 not estimated 
    % case 1                          
    opt.nutdy.int = opt.nutdy.int;         % case 1 interval in minutes  
    % case 1 
    opt.nutdy.constrain = opt.nutdy.constrain;      % case 1 if constraints should be applied
                                                    % 0 if not applied
        % case 1 
        opt.nutdy.coef = opt.nutdy.coef;      % constraint in mas/day 
                                              % (e.g. 0.0001 creates just one offset, 3 is loose)

% open opt-file
optname = ['../../../VLBI_OPT/',opt.diropt,'/',antenna(1).ngsfile(1:14),'.OPT'];
fid = fopen(optname);

for istat = 1:na
    opt.stat(istat).clkbreak = [];
    opt.stat(istat).delete = 0;
end

% if the OPT file is existing
if fid~=-1
    
    % check whether there is an entry for the clock reference in the OPT
    % file
    line = [fgetl(fid),'                        '];
    while (~strcmp(line(1:15),'CLOCK REFERENCE'))&&(~feof(fid))
        line = [fgetl(fid),'                        '];
    end
    if feof(fid)
        opt.stat(1).ref = 1;
    else
        line = [fgetl(fid) '        '];
        stat = line(1:8);
        for istat = 1:na
            if strcmp(antenna(istat).name,stat)
                opt.stat(istat).ref = 1;
            end
        end
    end
    
    frewind(fid);
    % check whether there is an entry for the clock breaks
    line = [fgetl(fid),'                        '];
    while (~strcmp(line(1:12),'CLOCK BREAKS'))&&(~feof(fid))
        line = [fgetl(fid),'                        '];
    end
    if feof(fid)
    else
        nbreak = sscanf(line(14:length(line)),'%i',[1,1]);  % Boehm 21 Aug 2009
        for j = 1:nbreak
            line = fgetl(fid);
            stat = line(1:8); 
            bmjd = sscanf(line(9:length(line)),'%f',[1,1]);  % Boehm 21 Aug 2009
            for istat = 1:na
                if strcmp(antenna(istat).name,stat)
                    opt.stat(istat).clkbreak = [opt.stat(istat).clkbreak; bmjd];
                end
            end
        end
    end
    
    if exist('stat','var') &  find([opt.stat.ref])==1 %25/06/2014
       if length(opt.stat) < 2
            if ~isempty(opt.stat(1).clkbreak) 
                opt.stat(1).ref = 0; opt.stat(2).ref = 1;
                if ~isempty(opt.stat(2).clkbreak)
                    opt.stat(1).ref = 0; opt.stat(3).ref = 1;
                end
            end
       end
    end

    
    
    frewind(fid);
    % check whether stations should be excluded
    line = [fgetl(fid),'                        '];
    while (~strcmp(line(1:23),'STATIONS TO BE EXCLUDED'))&&(~feof(fid))
        line = [fgetl(fid),'                        '];
    end
    if feof(fid)
    else
        ndel = str2num(line(25:length(line)));
        for j = 1:ndel
            line = fgetl(fid);
            stat = line(1:8); 
            for istat = 1:na
                if strcmp(antenna(istat).name,stat)
                    opt.stat(istat).delete = 1;
                end
            end
        end
    end
    
% if OPT file is not existing    
else
    opt.stat(1).ref = 1; % first clock of the session is seleced as reference clock
end

opt.ref_first_clk = find([opt.stat.ref]);

if fid~=-1
    fclose(fid);
end

% Outlier test options
opt.par_outlier = opt.par_outlier;  

opt.simple_outlier = opt.simple_outlier; % case 1 : par*mo is the boundary
                                         % value [fast]
 
                                         % case 0 : no outlier test

opt.basic_outlier = opt.basic_outlier;  % case 1 : par*mo*sqrt(qvv) is the boundary value [fast]
                                        % case 0 : no outlier test
                                        
% Estimation of parameters (hana 12Nov2010)
opt.est_singleses = opt.est_singleses;  % case 0: parameters will be not estimated, only N and b will be prepared
                                        
% GUI for vie_lsm
opt.control_gui_vie_lsm = opt.control_gui_vie_lsm;
if opt.control_gui_vie_lsm == 1
    vie_lsm_gui_first;
    uiwait(vie_lsm_gui_first);
end


