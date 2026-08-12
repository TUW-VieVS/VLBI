% ************************************************************************
%   Description:
%   function to find the reference epoch for the interpolation of the
%   satellite position
%
%   Reference:
%
%   Input:
%       scan_cur  - current scan
%       s_cur     - current satellite
%  
%
%   Output:
%       t_interp_normalized  - interpolation points
%       t_obs_normalized     - observations points
%       refidx               - reference index
%
%   External calls:
%
%   Coded for VieVS:
%   10 Aug 2026 by H. Wolf 
%
%   Revision:
%
% ************************************************************************

function [t_interp_normalized, t_obs_normalized, refidx] = get_reference_epoch_satellite(scan_cur, s_cur)

    continuous_time = s_cur.mjd * 86400;
    obs_time_continuous = scan_cur.mjd * 86400;
    [~, refidx_center] = min(abs(continuous_time - obs_time_continuous));
    start_idx = refidx_center - 7;
    end_idx = refidx_center + 6;

    if start_idx < 1 || end_idx > length(continuous_time)
        error('Satellite orbit data do not cover the required time! Add missing data to orbit data file!');
    end

    refidx = start_idx:end_idx;
    t_interp_points = continuous_time(refidx);

    t_offset = t_interp_points(1);
    t_interp_normalized = t_interp_points - t_offset;

    t_obs_normalized = obs_time_continuous - t_offset;
end