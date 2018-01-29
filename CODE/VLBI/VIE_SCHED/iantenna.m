% Purpose  
%   Read antenna information from the antenna catalog file.
% History  
%   2010-04-06   Jing SUN   Created
%
% CHANGES
% - 2015-08-24, A. Hellerschmied: The parameters PARA.MARGEL1 and PARA.MARGEL2 only influence the axis limts of AzEl- and HaDec-mounted antennas. 
% - 2016-06-27, M. Schartner: different computation to calc neutral range of cable wrap


function [station] = iantenna(filename, station, PARA)

% open antenna catalog file
fid = fopen(filename, 'r');
if (fid < 0)
    error('    no antenna catalog file %s !\n', filename);
end

% read the antenna information
stanum = length(station);
for ista = 1 : stanum
    frewind(fid);
    ifantenna = 0;
    while ~feof(fid)
        line = fgetl(fid);
        linelength = length(line);
        if ((linelength < 11) | (~strcmp(line(1), ' ')))   %%%%%
            continue;
        end
        if ~strcmp(line(4:11), station(ista).name(1:8))
            continue;
        end
        ifantenna = 1;
        station(ista).id   = line(2);
        station(ista).axis = line(13:16);
        [tmprd, count, errmsg, nextindex] = sscanf(line(18:linelength), '%f', 10);
        station(ista).offset = tmprd(1);
        station(ista).rate1  = tmprd(2) / 60;         % [deg/s]
        station(ista).c1     = tmprd(3);              % [sec] 
        if strcmpi(station(ista).axis, 'AZEL') || strcmpi(station(ista).axis, 'HADC') % Only HADC and AZEL mpunt types ==> Marge from param.txt is taken into account.
            station(ista).lim11  = (tmprd(4)+PARA.MARGEL1) * pi / 180.0; % [rad]
            station(ista).lim12  = (tmprd(5)-PARA.MARGEL1) * pi / 180.0; % [rad]
            station(ista).lim21  = (tmprd(8)+PARA.MARGEL2) * pi / 180.0; % [rad]
            station(ista).lim22  = (tmprd(9)-PARA.MARGEL2) * pi / 180.0; % [rad]
        else % Other antenna mount types (XYEW and XYNS, etc...) ==> Marge from param.txt is NOT taken into account.
            station(ista).lim11  = tmprd(4) * pi / 180.0; % [rad]
            station(ista).lim12  = tmprd(5) * pi / 180.0; % [rad]
            station(ista).lim21  = tmprd(8) * pi / 180.0; % [rad]
            station(ista).lim22  = tmprd(9) * pi / 180.0; % [rad]
        end
        station(ista).rate2  = tmprd(6) / 60;         % [deg/s]
        station(ista).c2     = tmprd(7);              % [sec] 
        station(ista).diam   = tmprd(10);             % [meter]
        index = 18 + nextindex - 1;
        [tmppo, count, errmsg, nextindex] = sscanf(line(index:linelength), '%s', 1);
        station(ista).po(1:2) = tmppo(1:2);
        index = index + nextindex - 1;
        [tmpeq, count, errmsg, nextindex] = sscanf(line(index:linelength), '%s', 1);
        station(ista).eq(1:length(tmpeq)) = tmpeq(1:length(tmpeq));
        index = index + nextindex - 1;
        [tmpms, count, errmsg, nextindex] = sscanf(line(index:linelength), '%s', 1);
        if (count == 1)
            station(ista).ms(1:2) = tmpms(1:2);
        end
        fprintf(PARA.fid_header,'    %s\n', line(1:linelength));
        break;
    end
    if (ifantenna == 0)
        error('    There is no antenna information for station %s \n', station(ista).name);
    end
    if (~strcmp(station(ista).axis(1:4), 'AZEL'))&(~strcmp(station(ista).axis(1:4), 'HADC'))&(~strcmp(station(ista).axis(1:4), 'XYEW'))
        error('    Axis type %s can not be processed\n', station(ista).axis(1:4));
    end
end

% close file
fclose(fid);

% cable wrap
for ista = 1 : stanum
    az1 = station(ista).lim11*180/pi - PARA.MARGEL1; %%% no marge
    az2 = station(ista).lim12*180/pi + PARA.MARGEL1; %%% no marge
    station(ista).aznp = (az2-az1)/2 + az1;
    % neutral range
    
    az = max([az1,az2-360]);
    station(ista).azn1 = az;
%     az = station(ista).aznp;
%     while ((az>=az1)&(az<=az2))
%         az = az - 0.001;
%         if ((az+360)<=az2) 
%             break;
%         end  
%     end
%     station(ista).azn1 = az + 0.001;
%     
%     az = station(ista).aznp;
%     while ((az>=az1)&(az<=az2))
%         az = az + 0.001;
%         if ((az-360)>=az1)
%             break;
%         end  
%     end
% station(ista).azn2 = az - 0.001;
    az = min([az1+360,az2]);
    station(ista).azn2 = az;
    if (station(ista).azn1 < station(ista).lim11*180/pi)
        station(ista).azn1 = station(ista).lim11*180/pi;
    end
    if (station(ista).azn2 > station(ista).lim12*180/pi)
        station(ista).azn2 = station(ista).lim12*180/pi;
    end
    % clockwise range
    station(ista).azc1 = station(ista).azn2;
    station(ista).azc2 = az2-PARA.MARGEL1;
    % count-clockwise range
    station(ista).azw1 = az1+PARA.MARGEL1;
    station(ista).azw2 = station(ista).azn1;
    % deg -> rad
    station(ista).aznp = station(ista).aznp * pi/180;
    station(ista).azn1 = station(ista).azn1 * pi/180;
    station(ista).azn2 = station(ista).azn2 * pi/180;
    station(ista).azc1 = station(ista).azc1 * pi/180;
    station(ista).azc2 = station(ista).azc2 * pi/180;
    station(ista).azw1 = station(ista).azw1 * pi/180;
    station(ista).azw2 = station(ista).azw2 * pi/180;
end


