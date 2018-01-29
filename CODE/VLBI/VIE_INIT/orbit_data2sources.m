% #########################################################################
% #     orbit_data2sources
% #########################################################################
%
% DESCRITPION
% This function adds the orbit information of the "orbit_data" structure to the 
% VieVS "sources" structure
%
% AUTHOR 
%   A. Hellerschmied, 2016-08-01
%
% INPUT
%  - orbit_data          - structure containing the orbit data
%  - sources             - VieVS sources structure
%  
% OUTPUT
%  - sources             - VieVS sources structure (updates satellite positions)
%
%
% CHANGES
%  - yyyy-mm-dd, <first + second name>: Description
%  - 2016-09-21, A. Hellerschmied: Added field vx_trf, vy_trf, vz_trf and flag_v_trf to orbit_data and source.s structures 
%  - 2016-12-05, A. Hellerschmied: Added fields 'sec_of_day', 'year', 'month', 'day', 'hour', 'minu', 'sec' to orbit_data and source.s structures 


function [sources] = orbit_data2sources(orbit_data, sources)

% loop over all satellites in sources:
for i_sat = 1 : length(sources.s)
    
    % Find matching orbit_data entry:
    orbit_data_ind = strcmp({orbit_data.sat.name}, sources.s(i_sat).name);
    sources.s(i_sat).x_trf      = orbit_data.sat(orbit_data_ind).x_trf;
    sources.s(i_sat).y_trf      = orbit_data.sat(orbit_data_ind).y_trf;
    sources.s(i_sat).z_trf      = orbit_data.sat(orbit_data_ind).z_trf;
    sources.s(i_sat).vx_trf     = orbit_data.sat(orbit_data_ind).vx_trf;
    sources.s(i_sat).vy_trf     = orbit_data.sat(orbit_data_ind).vy_trf;
    sources.s(i_sat).vz_trf     = orbit_data.sat(orbit_data_ind).vz_trf;
    sources.s(i_sat).flag_v_trf = orbit_data.sat(orbit_data_ind).flag_v_trf;
    sources.s(i_sat).mjd        = orbit_data.epoch_mjd;
    sources.s(i_sat).year       = orbit_data.year;
    sources.s(i_sat).month      = orbit_data.month;
    sources.s(i_sat).day        = orbit_data.day;
    sources.s(i_sat).hour       = orbit_data.hour;
    sources.s(i_sat).minu       = orbit_data.minu;
    sources.s(i_sat).sec        = orbit_data.sec;
    sources.s(i_sat).sec_of_day = orbit_data.sec_of_day;
    
end



