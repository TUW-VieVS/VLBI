function [Ax, Ay,Az, Ax_s, Ay_s, Az_s, Ax_q, Ay_q, Az_q] = a_xyz(obs_per_stat,nobserv,per_stat, opt)
    Ax(nobserv,1) = 0;
    Ay(nobserv,1) = 0;
    Az(nobserv,1) = 0;

    Ax_s(nobserv,1) = 0;
    Ay_s(nobserv,1) = 0;
    Az_s(nobserv,1) = 0;

    Ax_q(nobserv,1) = 0;
    Ay_q(nobserv,1) = 0;
    Az_q(nobserv,1) = 0;

    if opt.stc_sat == 1 || opt.stc_qs == 1
        if isfield(obs_per_stat, 'first_sat')
            first_s = obs_per_stat.first_sat; % -1 if it is i1 OR +1 if it is i2 
            dx_s = obs_per_stat.dx_sat; % partial derivatives of the observations with respect to dx (one specific station)
            dy_s = obs_per_stat.dy_sat; % partial derivatives of the observations with respect to dy (one specific station)
            dz_s = obs_per_stat.dz_sat; % partial derivatives of the observations with respect to dz (one specific station)
            nob_s = per_stat.oc_nob_sat;
        else
            error('ERROR: Session contains no Satellite Observations - Station Coordinate Estimation from Satellite Observations is therefore not possible!')
        end
        Ax_s(nob_s,1) = first_s.*dx_s;
        Ay_s(nob_s,1) = first_s.*dy_s;
        Az_s(nob_s,1) = first_s.*dz_s;
    else
        Ax_s = [];
        Ay_s = [];
        Az_s = [];
    end

    if opt.stc_qu == 1 || opt.stc_qs == 1
        if isfield(obs_per_stat, 'first_qu')
            first_q = obs_per_stat.first_qu; % -1 if it is i1 OR +1 if it is i2 
            dx_q = obs_per_stat.dx_qu; % partial derivatives of the observations with respect to dx (one specific station)
            dy_q = obs_per_stat.dy_qu; % partial derivatives of the observations with respect to dy (one specific station)
            dz_q = obs_per_stat.dz_qu; % partial derivatives of the observations with respect to dz (one specific station)
            nob_q = per_stat.oc_nob_qu;
        else
            error('ERROR: Session contains no Quasar Observations - Station Coordinate Estimation from Quasar Observations is therefore not possible!')
        end
        Ax_q(nob_q,1) = first_q.*dx_q;
        Ay_q(nob_q,1) = first_q.*dy_q;
        Az_q(nob_q,1) = first_q.*dz_q;
    else
        Ax_q = [];
        Ay_q = [];
        Az_q = [];
    end

    if opt.stc_all == 1
        first = obs_per_stat.first; % -1 if it is i1 OR +1 if it is i2 
        %total = obs_per_stat.total; % total number of observations that are carried out by the station
        %nob = obs_per_stat.nob; % the row numbers of the observations of that station in the oc vector, A, and P matrices 
        dx = obs_per_stat.dx; % partial derivatives of the observations with respect to dx (one specific station)
        dy = obs_per_stat.dy; % partial derivatives of the observations with respect to dy (one specific station)
        dz = obs_per_stat.dz; % partial derivatives of the observations with respect to dz (one specific station)
        nob = per_stat.oc_nob;

        % assigning the partial derivatives of station coordinates to the specific rows 
        Ax(nob,1) = first.*dx;
        Ay(nob,1) = first.*dy;
        Az(nob,1) = first.*dz; 
    else
        Ax = [];
        Ay = [];
        Az = [];
    end
end