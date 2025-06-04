% ************************************************************************
%   Description:
%	Organizing the data of satellite and storing it in a per_satellite struct
%
%
%   Input:										
%     opt                 opt settings
%     scan                scan struct
%     obs_per_satellite   observations per satellite
%     mjd0                beginning of the session
% 
% 
%   Output:
%     per_satellite       data (derivatives, intervals, ..) per satellite
% 
%   External calls: 	
%
%       
%   Coded for VieVS: 
%   17 Dec 2024 by Helene Wolf
%
%   Revision: 
%
%
% ************************************************************************

function [per_satellite] = get_per_satellite(opt, scan, obs_per_satellite, mjd0, name_orb)

    n_sat = length(opt.satellite);

    per_satellite =struct();

    for iSat = 1 : n_sat    
        mjd_sat = obs_per_satellite(iSat).mjd;
        per_satellite(iSat).name = opt.satellite(iSat).name;
        mjd1 = min(mjd_sat);    % time of the first scan of the current source in mjd [days]
        mjd2 = max(mjd_sat);    % time of the last scan of the current sourc in mjd [days]
        
        flag_within_interval = false;
        
        % mjd0(midnight)--------mjd1(start)------(Session)------mjd2(end)-----
        t1 = (mjd1-mjd0)*24*60; % time between first scan of this source and mjd0 [minutes] 
        t2 = (mjd2-mjd0)*24*60; % time between last scan of this source and mjd0 [minutes]
        
        if abs(t1 - round(t1)) < 0.0001
            t1 = round(t1);
        end
        if abs(t2 - round(t2)) < 0.0001
            t2 = round(t2);
        end
    
        if opt.SatPos.pw_sat
            % Estimation intervals of satellite pos.
            t10 = floor(t1/opt.SatPos.sat_pos_int) * opt.SatPos.sat_pos_int;
            t20 = ceil(t2/opt.SatPos.sat_pos_int)  * opt.SatPos.sat_pos_int;
       
            if opt.SatPos.sat_pos_int == 1440
                T = t1+(t2 - t1)/2; % one value
                per_satellite(iSat).est_int_id_sat_pos = 0;
            elseif (t2 - t1) < opt.SatPos.sat_pos_int +2 && (t2 - t1) > opt.SatPos.sat_pos_int -2
               T = [t1, t1+opt.SatPos.sat_pos_int];
               per_satellite(iSat).est_int_id_sat_pos = 0;
            else
                if t10 + opt.SatPos.sat_pos_int >= t20
                    t10 = t1;
                    t20 = t2;
                    T = [t10, t20];
                else
                    T = t10 : opt.SatPos.sat_pos_int : t20;   % Estimation epochs for sat. coor.
                end
            end             
            per_satellite(iSat).T_sat_pos = T;
            per_satellite(iSat).n_unk_sat_pos = length(T) - 1; % Number of estimation intervals for sat. coor.
        end


        if opt.KepEle.estKepEle  
            for iKepEle = 1:6
                if opt.KepEle.('estKepEle' + string(iKepEle))
                    % Estimation intervals of satellite pos.
                    int_min = opt.KepEle.('estIntKepEle' + string(iKepEle));
                    t10 = floor(t1/int_min) * int_min;
                    t20 = ceil(t2/int_min) * int_min;
                
                    if int_min == 1440
                        T = t1+(t2 - t1)/2; %1 value
                        per_satellite(iSat).('est_int_id_KepEle' + string(iKepEle)) = 0; 
                    elseif (t2 - t1) < int_min +2 && (t2 - t1) > int_min -2
                        T = [t1, t1+int_min];
                        per_satellite(iSat).('est_int_id_KepEle' + string(iKepEle))  = 0;
                    else
                        if t10 + int_min >= t20
                            t10 = t1;
                            t20 = t2;
                            T = [t10, t20];
                        else
                            T = t10 : int_min : t20;   % Estimation epochs for sat. coor.
                        end
                    end
                    per_satellite(iSat).('T_KepEle' + string(iKepEle)) = T;
                    per_satellite(iSat).('n_unk_KepEle' + string(iKepEle)) = length(T) - 1; % Number of estimation intervals for sat. coor.
                end
            end
        end

        i_obs_in_sess   = 0; % observation index in this session (absolut!)
        n_obs_of_sat    = 0; % observation index of the current source
        for iScan = 1 : opt.scans_total
            for iObs = 1 : scan(iScan).nobs    
                i_obs_in_sess = i_obs_in_sess + 1; 
                
                if strcmp(scan(iScan).obs_type, 's') 
                    if scan(iScan).iso == iSat
                        n_obs_of_sat = n_obs_of_sat + 1;
                        
                        per_satellite(iSat).mjd(n_obs_of_sat) = scan(iScan).mjd;     % time of the observations per source [day]
                        per_satellite(iSat).nob(n_obs_of_sat) = i_obs_in_sess;       % row number of observation in o-c per source
                        
                        if opt.SatPos.pw_sat
                            % Assign partial derivarives:
                            switch(opt.SatPos.sat_pos_est_ref_frame)
                                case 'gcrf'     % GCRF [sec/m]
                                    per_satellite(iSat).pd_pos1(n_obs_of_sat)    = scan(iScan).obs(iObs).psat_gcrf(1);
                                    per_satellite(iSat).pd_pos2(n_obs_of_sat)    = scan(iScan).obs(iObs).psat_gcrf(2);
                                    per_satellite(iSat).pd_pos3(n_obs_of_sat)    = scan(iScan).obs(iObs).psat_gcrf(3);
                                case 'trf'      % TRF [sec/m]
                                    per_satellite(iSat).pd_pos1(n_obs_of_sat)    = scan(iScan).obs(iObs).psat_trf(1);
                                    per_satellite(iSat).pd_pos2(n_obs_of_sat)    = scan(iScan).obs(iObs).psat_trf(2);
                                    per_satellite(iSat).pd_pos3(n_obs_of_sat)    = scan(iScan).obs(iObs).psat_trf(3);
                                case 'rsw'      % RSW system ("satellite coord. sys.") [sec/m]
                                    per_satellite(iSat).pd_pos1(n_obs_of_sat)    = scan(iScan).obs(iObs).psat_rsw(1);
                                    per_satellite(iSat).pd_pos2(n_obs_of_sat)    = scan(iScan).obs(iObs).psat_rsw(2);
                                    per_satellite(iSat).pd_pos3(n_obs_of_sat)    = scan(iScan).obs(iObs).psat_rsw(3);
                                case 'ntw'      % NTW system [sec/m]
                                    per_satellite(iSat).pd_pos1(n_obs_of_sat)    = scan(iScan).obs(iObs).psat_ntw(1);
                                    per_satellite(iSat).pd_pos2(n_obs_of_sat)    = scan(iScan).obs(iObs).psat_ntw(2);
                                    per_satellite(iSat).pd_pos3(n_obs_of_sat)    = scan(iScan).obs(iObs).psat_ntw(3);
                            end
                            per_satellite(iSat).minute_sat_pos(n_obs_of_sat) = (scan(iScan).mjd - mjd0) * 24*60; % time reference in minutes since epoch mjd0
                            if abs(per_satellite(iSat).minute_sat_pos(n_obs_of_sat) - round(per_satellite(iSat).minute_sat_pos(n_obs_of_sat))) < 0.0001
                                per_satellite(iSat).minute_sat_pos(n_obs_of_sat)  = round(per_satellite(iSat).minute_sat_pos(n_obs_of_sat) );
                            end

                            % Loop over all estimation intervals of the current obs. target:
                            for i_int = 1 : per_satellite(iSat).n_unk_sat_pos
                                if i_int == per_satellite(iSat).n_unk_sat_pos  % Last interval
                                    flag_within_interval = (per_satellite(iSat).minute_sat_pos(n_obs_of_sat) >= per_satellite(iSat).T_sat_pos(i_int)) && (per_satellite(iSat).minute_sat_pos(n_obs_of_sat) <= per_satellite(iSat).T_sat_pos(i_int+1));
                                else                % All other intervals
                                    flag_within_interval = (per_satellite(iSat).minute_sat_pos(n_obs_of_sat) >= per_satellite(iSat).T_sat_pos(i_int)) && (per_satellite(iSat).minute_sat_pos(n_obs_of_sat) < per_satellite(iSat).T_sat_pos(i_int+1));
                                end
                                if flag_within_interval
                                    per_satellite(iSat).est_int_id_sat_pos(n_obs_of_sat) = i_int; % Index of estimation interval that contains the current observation of the source
                                end
                            end

                       end

                       if opt.KepEle.estKepEle
                            for iKepEle = 1:6
                                if opt.KepEle.('estKepEle' + string(iKepEle))
                                    per_satellite(iSat).('pd_KepEle' + string(iKepEle))(n_obs_of_sat) = scan(iScan).obs(iObs).(name_orb)(iKepEle);
                                
                                    per_satellite(iSat).('minute_KepEle' + string(iKepEle))(n_obs_of_sat) = (scan(iScan).mjd - mjd0) * 24*60; % time reference in minutes since epoch mjd0
                                    if abs(per_satellite(iSat).('minute_KepEle' + string(iKepEle))(n_obs_of_sat) - round(per_satellite(iSat).('minute_KepEle' + string(iKepEle))(n_obs_of_sat))) < 0.0001
                                        per_satellite(iSat).('minute_KepEle' + string(iKepEle))(n_obs_of_sat)  = round(per_satellite(iSat).('minute_KepEle' + string(iKepEle))(n_obs_of_sat) );
                                    end
                            
                                    % Loop over all estimation intervals of the current obs. target:
                                    for i_int = 1 : per_satellite(iSat).('n_unk_KepEle' + string(iKepEle))
                                        if i_int == per_satellite(iSat).('n_unk_KepEle' + string(iKepEle))  % Last interval
                                            flag_within_interval = (per_satellite(iSat).('minute_KepEle' + string(iKepEle))(n_obs_of_sat) >= per_satellite(iSat).('T_KepEle' + string(iKepEle))(i_int)) && (per_satellite(iSat).('minute_KepEle' + string(iKepEle))(n_obs_of_sat) <= per_satellite(iSat).('T_KepEle' + string(iKepEle))(i_int+1));
                                        else                % All other intervals
                                            flag_within_interval = (per_satellite(iSat).('minute_KepEle' + string(iKepEle))(n_obs_of_sat) >= per_satellite(iSat).('T_KepEle' + string(iKepEle))(i_int)) && (per_satellite(iSat).('minute_KepEle' + string(iKepEle))(n_obs_of_sat) < per_satellite(iSat).('T_KepEle' + string(iKepEle))(i_int+1));
                                        end
                                        if flag_within_interval
                                            per_satellite(iSat).('est_int_id_KepEle' + string(iKepEle))(n_obs_of_sat) = i_int; % Index of estimation interval that contains the current observation of the source
                                        end
                                    end
                                end
                            end
                        end
                    end 
                end 
            end 
        end     
        per_satellite(iSat).total = n_obs_of_sat;
    end
end