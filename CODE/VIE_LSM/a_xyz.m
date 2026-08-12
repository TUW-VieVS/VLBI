function [Ax, Ay,Az, Ax_s, Ay_s, Az_s, Ax_q, Ay_q, Az_q] = a_xyz(obs_per_stat,nobserv,per_stat, opt)

    Ax = zeros(nobserv, 1); Ay = zeros(nobserv, 1); Az = zeros(nobserv, 1);
    Ax_s = zeros(nobserv, 1); Ay_s = zeros(nobserv, 1); Az_s = zeros(nobserv, 1);
    Ax_q = zeros(nobserv, 1); Ay_q = zeros(nobserv, 1); Az_q = zeros(nobserv, 1);

    if opt.stc_sat == 1 || opt.stc_qs == 1
        if isfield(obs_per_stat, 'first_sat')
            nob_s = per_stat.oc_nob_sat;
            Ax_s(nob_s,1) = obs_per_stat.first_sat .* obs_per_stat.dx_sat;
            Ay_s(nob_s,1) = obs_per_stat.first_sat .* obs_per_stat.dy_sat;
            Az_s(nob_s,1) = obs_per_stat.first_sat .* obs_per_stat.dz_sat;
        else
            error('ERROR: Session contains no Satellite Observations! - Station Coordinate Estimation from Satellite Observations is therefore not possible!')
        end
    else 
        Ax_s = []; Ay_s = []; Az_s = [];
    end

    if opt.stc_qu == 1 || opt.stc_qs == 1
        if isfield(obs_per_stat, 'first_qu')
            nob_q = per_stat.oc_nob_qu;
            Ax_q(nob_q,1) = obs_per_stat.first_qu .* obs_per_stat.dx_qu;
            Ay_q(nob_q,1) = obs_per_stat.first_qu .* obs_per_stat.dy_qu;
            Az_q(nob_q,1) = obs_per_stat.first_qu .* obs_per_stat.dz_qu;
        else
            error('ERROR: Session contains no Quasar Observations! - Station Coordinate Estimation from Quasar Observations is therefore not possible!')
        end
    else 
        Ax_q = []; Ay_q = []; Az_q = [];
    end

    if opt.stc_all == 1
         nob = per_stat.oc_nob;
        Ax(nob,1) = obs_per_stat.first .* obs_per_stat.dx;
        Ay(nob,1) = obs_per_stat.first .* obs_per_stat.dy;
        Az(nob,1) = obs_per_stat.first .* obs_per_stat.dz;
    else
        Ax = []; Ay = []; Az = [];
    end
end