% ************************************************************************
%   Description:
%   function to rearrange the data from scan-based to satellite-based
%
%   Reference:
%
%   Input:
%       'opt'                   structure array             (for info. /DOC/opt.doc)
%       'scan'                  structure array             (for info. /DOC/scan.doc)
%       'num_of_scans'          (1,1)                       number of scans in a session
%       'mjd_sat'               (num. of obs. to source,1)  epochs of the observations to the source in mjd
%       'i_sat'                 (1,1)                       Satellite index (ref. to sources.s)
%       'mjd0'                  (1,1)                       mjd of the 0:00 UTC of the first day of the session
%       'sat_pos_int_minutes'   (1,optional)                estimation intervals of the satellite coordinates
%
%   Output:
%       'per_satellite'    structure array     data for the currents satellite (with index i_sat)
%       'n_unk'            structure array     number of estimation intervals
%       'T_'               structure array     estimation epochs
%
%   External calls:
%
%   Coded for VieVS:
%   03 May 2017 by A. Hellerschmied
%
%   Revisions:
%
% ************************************************************************
function [per_satellite, n_unk, T_] = satellitewisepar(opt, scan, num_of_scans, mjd_sat, i_sat, mjd0, sat_pos_int_minutes)

mjd1 = min(mjd_sat);    % time of the first scan of the current source in mjd [days]
mjd2 = max(mjd_sat);    % time of the last scan of the current sourc in mjd [days]

flag_within_interval = false;

% mjd0(midnight)--------mjd1(start)------(Session)------mjd2(end)-----
t1 = (mjd1-mjd0)*24*60; % time between first scan of this source and mjd0 [minutes] 
t2 = (mjd2-mjd0)*24*60; % time between last scan of this source and mjd0 [minutes] 

% Estimation intervals of satellite pos.
t10 = floor(t1/sat_pos_int_minutes) * sat_pos_int_minutes;
t20 = ceil(t2/sat_pos_int_minutes) * sat_pos_int_minutes;

T_.sat = t10 : sat_pos_int_minutes : t20;       % Estimation epochs for sat. coor.
n_unk.sat = length(T_.sat) - 1;                 % Number of estimation intervals for sat. coor.

per_satellite.name = opt.satellite(i_sat).name; % The name of the source

% loop init.:
i_obs_in_sess   = 0; % observation index in this session (absolut!)
n_obs_of_sat    = 0; % observation index of the current source

% Loop over all scans:
for i_scan = 1 : num_of_scans
    
    % Loop over all observations in scan:
    for i_obs = 1 : scan(i_scan).nobs
        
        i_obs_in_sess = i_obs_in_sess + 1; 
        
        % Only consider observations of satellites:
        if strcmp(scan(i_scan).obs_type, 's') 
            
            if scan(i_scan).iso == i_sat
                
                n_obs_of_sat = n_obs_of_sat + 1;
                
                per_satellite.mjd(n_obs_of_sat) = scan(i_scan).mjd;     % time of the observations per source [day]
                per_satellite.nob(n_obs_of_sat) = i_obs_in_sess;        % row number of observation in o-c per source
                
                % Assign partial derivarives:
                switch(opt.sat_pos_est_ref_frame)
                    case 'gcrf'     % GCRF [sec/m]
                        per_satellite.pd_pos1(n_obs_of_sat)    = scan(i_scan).obs(i_obs).psou_sat_gcrf(1);
                        per_satellite.pd_pos2(n_obs_of_sat)    = scan(i_scan).obs(i_obs).psou_sat_gcrf(2);
                        per_satellite.pd_pos3(n_obs_of_sat)    = scan(i_scan).obs(i_obs).psou_sat_gcrf(3);
                    case 'trf'      % TRF [sec/m]
                        per_satellite.pd_pos1(n_obs_of_sat)    = scan(i_scan).obs(i_obs).psou_sat_trf(1);
                        per_satellite.pd_pos2(n_obs_of_sat)    = scan(i_scan).obs(i_obs).psou_sat_trf(2);
                        per_satellite.pd_pos3(n_obs_of_sat)    = scan(i_scan).obs(i_obs).psou_sat_trf(3);
                    case 'rsw'      % RSW system ("satellite coord. sys.") [sec/m]
                        per_satellite.pd_pos1(n_obs_of_sat)    = scan(i_scan).obs(i_obs).psou_sat_rsw(1);
                        per_satellite.pd_pos2(n_obs_of_sat)    = scan(i_scan).obs(i_obs).psou_sat_rsw(2);
                        per_satellite.pd_pos3(n_obs_of_sat)    = scan(i_scan).obs(i_obs).psou_sat_rsw(3);
                end
        
                per_satellite.minute(n_obs_of_sat) = (scan(i_scan).mjd - mjd0) * 24*60; % time reference in minutes since epoch mjd0
                
                % Loop over all estimation intervals of the current obs. target:
                for i_int = 1 : n_unk.sat
                    if i_int == n_unk.sat  % Last interval
                        flag_within_interval = (per_satellite.minute(n_obs_of_sat) >= T_.sat(i_int)) && (per_satellite.minute(n_obs_of_sat) <= T_.sat(i_int+1));
                    else                % All other intervals
                        flag_within_interval = (per_satellite.minute(n_obs_of_sat) >= T_.sat(i_int)) && (per_satellite.minute(n_obs_of_sat) < T_.sat(i_int+1));
                    end
                    if flag_within_interval
                        per_satellite.est_int_id_sat_pos(n_obs_of_sat) = i_int; % Index of estimation interval that contains the current observation of the source
                    end
                end
            end % if scan(i_scan).iso == i_sat
        end % if strcmp(scan(i_scan).obs_type, 's')
    end % for iobs = 1:scan(i_scan).nobs
end % for i_scan = 1 : num_of_scans

per_satellite.total = n_obs_of_sat;
