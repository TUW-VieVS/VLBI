% ************************************************************************
%   Description:
%	Manipulates the position of the satellite by changing one of the 
%   orbital elements (extending semi-major axis, increasing inclination) 
%   for orbit data from a TLE file.
%
%   Input:										
%     sources               sources struct
%     parameter             parameter struct
%     T2C_s                 transformation matrix to convert from TRS to CRS
%     numKepEle             number of keplerian element
%							 1: semi-major axis
%                            2: eccentricity
%                            3: inclination
%                            4: argument of perigee
%                            5: right ascension of ascending node
%                            6: argument of latitude
%
% 
%   Output:
%     sourcesChanged       changed sources struct
%     dKepEle			   change of orbital element
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

function [sourcesChanged_KepEle, dKepEle] = manipulateTLEData(sources, parameter, T2C_s, numKepEle)

    
    sources_c = sources;
    mjd_firstSatObs = min([sources.firstObsMjd]);
    mjd_lastSatObs = max([sources.lastObsMjd]);
    jd_firstSatObs = mjd_firstSatObs +  2400000.5;
    jd_lastSatObs = mjd_lastSatObs +  2400000.5;
    str_ind = max([strfind(parameter.vie_init.sc_orbit_file_path_name, '\'), strfind(parameter.vie_init.sc_orbit_file_path_name, '/')]);
    sat_orbit_file_path = parameter.vie_init.sc_orbit_file_path_name(1 : str_ind); 
    sat_orbit_file_name = parameter.vie_init.sc_orbit_file_path_name(str_ind+1 : end);
    [TLE_C, ~, ~ ] = read_tle(sat_orbit_file_path, sat_orbit_file_name);

    switch numKepEle
        case 1 %semimajor axis
            dKepEle = 1; %m 
            strTLE = 53; 
            strTLE2 = 64;
            isAngle = 0;
            f = '%011.8f';
        case 2 %eccentricity
            dKepEle = 0.000777;% []
            strTLE = 27; 
            strTLE2 = 34;
            isAngle = 0;
            f = '%09.7f';
        case 3 %inclination
            dKepEle = 1; %degrees
            strTLE = 9; 
            strTLE2 = 17;
            isAngle = 0;
            f = '%08.4f';
        case 4 %RAAN
            dKepEle = 1; %degrees
            strTLE = 18; 
            strTLE2 = 26; 
            isAngle = 1;
            f = '%08.4f';
        case 5 %argument of perigee (omega)
            dKepEle = 1; %degrees
            strTLE = 35; 
            strTLE2 = 43;
            isAngle = 1;
            f = '%08.4f'; 
        case 6 %argument of latitude
            dKepEle = 1; %degrees
            strTLE = 44; 
            strTLE2 = 52;
            isAngle = 1;
            f = '%08.4f';
    end

    TLE_KepEle = TLE_C;
    interval = 5; % [min]
    for i=1:length(TLE_C.tle_data.sat)
        line2 = TLE_C.tle_data.sat(i).line_2_str;
        if numKepEle ==1 
            mu = 3.9860044188 * 10^14; %m3s−2
            n = str2double(line2(53:64));
            a_o = mu^(1/3) / ((2*n*pi)/(86400))^(2/3);
            a_n = a_o + dKepEle;
            KepEle_n = ((1/a_n * mu^(1/3))^(3/2) * 86400) / (2*pi);

        elseif numKepEle == 2
           KepEle_o = str2double(append('0.', line2(strTLE:strTLE2)));
           KepEle_n = KepEle_o + dKepEle;
           test = num2str(KepEle_n', '%09.7f');
           KepEle_n = test(3:end);
           
        elseif numKepEle == 5
            KepEle_o = str2double(line2(strTLE:strTLE2));
            KepEle_n = KepEle_o + dKepEle; % rotate omega forward
            if KepEle_n<0
                KepEle_n = KepEle_n+360;
            end
        
            M = str2double(line2(44:52));
            e = str2double(append('0.', line2(27:34))); %eccentricity
            eps0=M; 
            eps1 = M + e*sind(eps0);
            while abs(eps0-eps1)>1e-10
                eps = M + e*sind(eps1);
                eps0=eps1;
                eps1= eps;
            end
            v = 2*atand(sqrt((1+e)/(1-e))*tand(eps/2));

            v_n = v - dKepEle; % rotate true anomaly backwards
            eps_n = 2 * atand( tand(v_n/2) * sqrt((1-e)/(1+e)) );
            M_n = eps_n-e*sind(eps_n);

            if M_n<0
                M_n = M_n +360;
            end

            TLE_KepEle.tle_data.sat(i).line_2_str = append(line2(1:44-1) , num2str(M_n, '%08.4f') , line2(52:end));
            line2 = TLE_KepEle.tle_data.sat(i).line_2_str;
        elseif numKepEle == 6
             M = str2double(line2(strTLE:strTLE2));
             e = str2double(append('0.', line2(27:34))); %eccentricity
             eps0=M; 
             eps1 = M + e*sind(eps0);
             while abs(eps0-eps1)>1e-10
                eps = M + e*sind(eps1);
                eps0=eps1;
                eps1= eps;
             end
             v = 2*atand(sqrt((1+e)/(1-e))*tand(eps/2));
             v_n = v+dKepEle;

             eps_n = 2 * atand( tand(v_n/2) * sqrt((1-e)/(1+e)) );
             M_n = eps_n-e*sind(eps_n);
             if M_n<0
                 M_n=M_n+360;
             end
             KepEle_n = M_n;
         else
           KepEle_o = str2double(line2(strTLE:strTLE2));
           KepEle_n = KepEle_o + dKepEle; 
        end

        if isAngle 
            if KepEle_n>360
                KepEle_n = KepEle_n - 360;
            end
            KepEle_n = num2str(KepEle_n', '%08.4f');         
        end
        TLE_KepEle.tle_data.sat(i).line_2_str = append(line2(1:strTLE-1) , num2str(KepEle_n, f) , line2(strTLE2:end));
    end

    sources_ch.s(1) = sources_c;
    [sat_data_KepEle, ~, ~] = tle_propagation(jd_firstSatObs - 2/24, jd_lastSatObs + 2/24, interval, TLE_KepEle);
    [sources_KepEle] = sat_data2sources(sat_data_KepEle, sources_ch); % position in crf, velocity in crf
    [sourcesChanged_KepEle] = getTRF_PosVelSatellites(sources_KepEle, T2C_s);

    if numKepEle>=3
        dKepEle = deg2rad(dKepEle); % degrees to rad
    end

    sourcesChanged_KepEle = sourcesChanged_KepEle.s(1);
end