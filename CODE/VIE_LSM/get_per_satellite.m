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
%   28-05-2026: changed the start and end to the session start and end
%   rather than first and last satellite observation; code optimization
%
% ************************************************************************

function [per_satellite] = get_per_satellite(opt, scan, mjd0, name_orb)
    obs_types      = {scan.obs_type};
    src_indices    = [scan.iso];
    is_sat         = strcmp(obs_types, 's');
    nsat = length(opt.satellite);
    per_satellite =struct();

    n_scans_total = length(scan);
    scan_start_global = zeros(n_scans_total, 1);
    if n_scans_total > 0
        scan_start_global(1) = 1;
        for k = 2:n_scans_total
            scan_start_global(k) = scan_start_global(k-1) + scan(k-1).nobs;
        end
    end

    mjd_start = scan(1).mjd;
    mjd_end = scan(end).mjd;
    session_duration = ceil((mjd_end - mjd_start)*24*60);

    est_mask_orb = [opt.ORB.params(:).estimate];
    est_mask_srp = [opt.SRP.params(:).estimate];
    est_orb  = find(est_mask_orb);
    est_srp  = find(est_mask_srp);

    tstart = (mjd_start - mjd0)*24*60;
    tend = ceil((mjd_end - mjd0)*24*60);
    

    for isat = 1 : nsat
        is_target_src  = ismember(src_indices, isat);
        mask = is_sat & is_target_src;
        orig_scan_idx = find(mask);
        selected_scans = scan(mask);

        per_satellite(isat).name = opt.satellite(isat).name;
                  
        if opt.SatPos.pw_sat
            int_min = opt.SatPos.sat_pos_int;
       
            if int_min == session_duration
                T = (tend-tstart)/2; %one value at mid of session
            elseif tstart + int_min >= tend
                T = [tstart, tend];
            else
                T = tstart : int_min : tend;  
            end             
            per_satellite(isat).T_pos = T;
            per_satellite(isat).n_unk_pos = length(T); % number of piecewise linear offsets
        end

        if opt.ORB.estORB
            for iorb = est_orb
                int_min = opt.ORB.params(iorb).interval; 
                
                if int_min == session_duration || int_min == 0 
                    T = (tend-tstart)/2; %one value at mid of session
                elseif tstart + int_min >= tend
                    T = [tstart, tend];
                else
                    T = tstart : int_min : tend;  
                end
                per_satellite(isat).('T_orb' + string(iorb)) = T;
                per_satellite(isat).('n_unk_orb' + string(iorb)) = length(T); % number of offsets
            end
        end

        if opt.SRP.estSRP  
            for isrp = est_srp
                int_min = opt.SRP.params(isrp).interval;

                if int_min == session_duration || int_min == 0 
                   T = (tend-tstart)/2; %one value at mid of session 
                elseif tstart + int_min >= tend
                    T = [tstart, tend];
                else    
                   T = tstart : int_min : tend; 
                end
                per_satellite(isat).('T_srp' + string(isrp)) = T;
                per_satellite(isat).('n_unk_srp' + string(isrp)) = length(T); % number of offsets
            end
        end

        nobs_sat    = 0; % observation index of the current source
         
        for k = 1:numel(selected_scans)
            s_cur = selected_scans(k);
            base_global_idx = scan_start_global(orig_scan_idx(k));

            for iObs = 1 : s_cur.nobs    
                nobs_sat = nobs_sat + 1;
                
                per_satellite(isat).mjd(nobs_sat) = s_cur.mjd;     % time of the observations per source [day]
                per_satellite(isat).nob(nobs_sat) = base_global_idx + iObs - 1;
                      
                t_minutes   = (s_cur.mjd - mjd0) * 24 * 60;
                round_mask = abs(t_minutes - round(t_minutes)) < 1e-4;
                t_minutes(round_mask) = round(t_minutes(round_mask));
                per_satellite(isat).minute(nobs_sat) = t_minutes;

                if opt.SatPos.pw_sat
                    
                    T = per_satellite(isat).T_pos;
                    n_unk = numel(T);
                    switch opt.SatPos.sat_pos_est_ref_frame
                        case 'gcrf', obs_fld = 'psat_gcrf';
                        case 'trf',  obs_fld = 'psat_trf';
                        case 'rsw',  obs_fld = 'psat_rsw';
                        case 'ntw',  obs_fld = 'psat_ntw';
                        otherwise,  error('Undefined Reference Frame: %s', opt.SatPos.sat_pos_est_ref_frame);
                    end
                    pos_vec = s_cur.obs(iObs).(obs_fld);
                    per_satellite(isat).pd_pos1(nobs_sat) = pos_vec(1);
                    per_satellite(isat).pd_pos2(nobs_sat) = pos_vec(2);
                    per_satellite(isat).pd_pos3(nobs_sat) = pos_vec(3);
                    
                    if n_unk > 1
                        int_ids = discretize(t_minutes, T);
                    else
                        int_ids = 1;
                    end
                    per_satellite(isat).est_int_id_pos(nobs_sat) = int_ids;
                end  
   
                if opt.ORB.estORB
                    for iorb = est_orb
                        T = per_satellite(isat).('T_orb' + string(iorb));
                        n_unk = numel(T);
                        raw_pd      = s_cur.obs(iObs).(name_orb)(iorb);
                        per_satellite(isat).('pd_orb' + string(iorb))(nobs_sat) = raw_pd;
                       
                        if n_unk > 1
                            int_ids = discretize(t_minutes, T);
                        else
                            int_ids = 1;
                        end
                        per_satellite(isat).('est_int_id_orb'+ string(iorb))(nobs_sat) = int_ids;
                    end
                 end

                if opt.SRP.estSRP
                    for isrp = est_srp
                        T = per_satellite(isat).('T_srp' + string(isrp));
                        n_unk = numel(T);
                        raw_pd      = s_cur.obs(iObs).('dsrp')(isrp); 
                        
                        per_satellite(isat).('pd_srp' + string(isrp))(nobs_sat) = raw_pd;
     
                        if n_unk > 1
                            int_ids = discretize(t_minutes, T);
                        else
                            int_ids = 1;
                        end
                        per_satellite(isat).('est_int_id_srp'+ string(isrp))(nobs_sat) = int_ids;           
                    end
                end
            end
        end 
        per_satellite(isat).total = nobs_sat;
    end
end