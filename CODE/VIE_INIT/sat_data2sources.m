% ************************************************************************
%   Description:
%   This function adds the orbit information from the "sat_data" structure to the 
%   VieVS "sources" structure
%   Input:										
%      sat_data          - structure containing the orbit data
%      sources           - VieVS sources structure
% 
%   Output:
%      sources           - VieVS sources structure (updates satellite positions)
% 
%   External calls: 	
%       
%       
%   Coded for VieVS: 
%   28 Apr 2022 by Helene Wolf
%
%   Revision: 
%  - 2024-12-13, Helene Wolf: changed the sat name to use the full name
%
% ************************************************************************

function [sources] = sat_data2sources(sat_data, sources)

    sat_names = strings(length(sat_data.sat),1);
    ids = strings(length(sat_data.sat),1);
    for i=1:length(sat_data.sat)
        sat_names(i) = regexprep(string(sat_data.sat(i).TLE_header_line(1:end)), ' ', '_');
        ids(i) = string(sat_data.sat(i).sat_number);
    end
    
    for iSat = 1 : length(sources.s)
        id = char(sources.s(iSat).id);
        id = id(1:end-1);
        orbit_data_ind = find(strcmp(ids, string(id)) == 1);
        if isempty(orbit_data_ind)
            error(fprintf(' *** Satellite %s is not included in TLE file.', regexprep(string(sources.s(iSat).name), '_', ' ') ))
        end
        r = [sat_data.sat(orbit_data_ind).prop(:).r];
        r = reshape(r,3,[]);
        r = r';
        v = [sat_data.sat(orbit_data_ind).prop(:).v];
        v = reshape(v,3,[]);
        v = v';
        sources.s(iSat).x_crf      = r(:, 1)*1000;
        sources.s(iSat).y_crf      = r(:, 2)*1000;
        sources.s(iSat).z_crf      = r(:, 3)*1000;
        sources.s(iSat).vx_crf     = v(:, 1)*1000;
        sources.s(iSat).vy_crf     = v(:, 2)*1000;
        sources.s(iSat).vz_crf     = v(:, 3)*1000;
        sources.s(iSat).flag_v_crf = 1;
        sources.s(iSat).mjd        = [sat_data.sat(orbit_data_ind).prop.jd]' - 2400000.5 .* ones(length([sat_data.sat(orbit_data_ind).prop.jd]),1);
        sources.s(iSat).year       = [sat_data.sat(orbit_data_ind).prop.year]';
        sources.s(iSat).month      = [sat_data.sat(orbit_data_ind).prop.mon]';
        sources.s(iSat).day        = [sat_data.sat(orbit_data_ind).prop.day]';
        sources.s(iSat).hour       = [sat_data.sat(orbit_data_ind).prop.h]';
        sources.s(iSat).minu       = [sat_data.sat(orbit_data_ind).prop.min]';
        sources.s(iSat).sec        = [sat_data.sat(orbit_data_ind).prop.sec]';
        sources.s(iSat).sec_of_day = [sat_data.sat(orbit_data_ind).prop.h]' .*3600 + [sat_data.sat(orbit_data_ind).prop.min]' .* 60 + [sat_data.sat(orbit_data_ind).prop.sec]';
        sources.s(iSat).tosc = fix(sources.s(iSat).firstObsMjd);
    end
end