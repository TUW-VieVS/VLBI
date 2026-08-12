% ************************************************************************
%   Description:
%	Computes the position and velocity for satellites based on the fso data
%   in a 5 minute interval. The position and velocity is then stored in the 
%   sources.s data for each satellite included in the NGS file.
%	The calculation is done using GPS Time - but in the sources.s struct the
%   the corresponding UTC time is stored.
%
%   Input:										
%     fso_data              fso_data
%     sources               sources struct
%
% 
%   Output:
%     sources              sources struct
% 
%   External calls: 	
%       getorb_par.m   datetime.m   modjuldat.m   mjd2date.m
%       
%   Coded for VieVS: 
%   17 Dec 2024 by Helene Wolf
%
%   Revision: 
%
%
% ************************************************************************

function [sources] = fso_data2sources(fso_data,sources)  
    ns = length(sources.s);

    for i_sat=1:ns
        prn = sources.s(i_sat).fso_name;

        interval = 5; % [min]
        time_gps = fso_data.tbound(1) : interval/1440 :fso_data.tbound(end); % in 5 minute steps
        [pos,vel] = getorb_par(prn,time_gps,fso_data.coeff,fso_data.t0,fso_data.hstep,fso_data.tbound,fso_data.satnum);
    
        dt = datetime(time_gps,'Format','yyyy MMMM dd hh:mm:ss.SSSSSSS','ConvertFrom','MJD');
        [year, month, day, hour, minute, sec] = datevec(dt);
    
        % Convert GPS time to UTC: get time diff. between UTC and GPS time:
        % tgps = UTC + leap_sec_tai_utc - 19 sec  => UTC = tgps + 19 sec - leap_sec_tai_utc
        leap_sec_tai_utc = tai_utc(time_gps);
        leap_sec_gps_utc = 19 - leap_sec_tai_utc;
        sec = sec + leap_sec_gps_utc';
            
        % Convert epochs to MJD:
        time_utc_again = modjuldat(year, month, day, hour, minute, sec);
        [year, month, day, hour, minute, sec] = mjd2date(time_utc_again); % round to integer seconds!
     
        sources.s(i_sat).x_crf = pos(:,1);
        sources.s(i_sat).y_crf = pos(:,2);
        sources.s(i_sat).z_crf = pos(:,3);
        sources.s(i_sat).vx_crf = vel(:,1);
        sources.s(i_sat).vy_crf = vel(:,2);
        sources.s(i_sat).vz_crf = vel(:,3);
        
        sources.s(i_sat).year = year';
        sources.s(i_sat).month = month';
        sources.s(i_sat).day = day';
        sources.s(i_sat).hour = hour';
        sources.s(i_sat).minu = minute';
        sources.s(i_sat).sec = sec';
        sources.s(i_sat).mjd = time_utc_again';
        sources.s(i_sat).mjd_gps = time_gps';
    
        sources.s(i_sat).sec_of_day = mod(time_utc_again, 1)'*24*60*60;
        sources.s(i_sat).flag_v_crf = 1;
        sources.s(i_sat).flag_v_trf = 0;

        sources.s(i_sat).tosc = fso_data.tosc;
    end
end