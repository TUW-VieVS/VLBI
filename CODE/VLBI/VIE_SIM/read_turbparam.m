% ************************************************************************
%   Description:
%   This function reads the parameters needed for the simulation from a
%   parameter file and stores them to a structure array.
%
%   References: 
%   ---
%
%   Input:		
%      fname:      string containing the name of the parameter file
%      antenna:    the VieVS antenna structure array (LEVEL1 data)
% 
%   Output:
%      a structure array containing the parameters for the simulation
%      dimension: [1 x # of antennas], the stations are in the same order
%      as in the antenna structure array!
% 
%   External calls: 
%   ---
%       
%   Coded for VieVS: 
%   July 2010 by Andrea Pany
%
%   Revision: 
%   2016-10-19, A. Hellerschmied: - Added the possibility to read wn separately for satellite observations (simparam(i).wn_sat).
%                                 - Availability of unambiguous entries for all stations in loaded turb file is now checked
% ************************************************************************

function simparam = read_turbparam(fname,antenna)

% Init.:
flag_wn_sat_available   = false;

% opfen parameter file
fid = fopen(fname);

% pre-allocate output
simparam = struct;

% skip header
fgetl(fid);

% read data from file
while ~feof(fid)
    line = fgetl(fid);
    for i = 1 : length(antenna)
        if strcmp(strtrim(line(1:8)),strtrim(antenna(i).name))
            vec   = sscanf(line(9:length(line)),'%f');
            simparam(i).Cn    = vec(1);     % [1e-7*m^-(1/3)]   refractive index structure constant
            simparam(i).H     = vec(2);     % [m]               effective height
            simparam(i).vn    = vec(3);     % [m/s]             north component of wind vector
            simparam(i).ve    = vec(4);     % [m/s]             east component of wind vector
            simparam(i).wzd0  = vec(5);     % [mm]              initial zenith wet delay
            simparam(i).dhseg = vec(6);     % [h]               time interval over which zwds are going to be correlated
            simparam(i).dh    = vec(7);     % [m]               height increment for integration
            simparam(i).sy1   = vec(8);     % []                ASD sy1 @ sy2 min
            simparam(i).sy2   = vec(9);     % [min]             ASD sy1 @ sy2 min
            simparam(i).wn    = vec(10);    % [ps]              white noise
            if length(vec) > 10
                simparam(i).wn_sat    = vec(11);    % [ps]      white noise for satellite observations
                flag_wn_sat_available = true;
            end
        end
    end
end

% Check, if one entry in the turb file was found for each station:
if length(simparam) < length(antenna)
    error('No entries for at least one station in the turb file %s!\n', fname);
elseif length(simparam) > length(antenna)
    error('Multiple entries for a single station in the turb file %s!\n', fname);
end

if flag_wn_sat_available
    fprintf('White noise for satellite observations defined in turbulance file: %s\n', fname)
    for i = 1 : length(antenna)
        fprintf('  - %8s : %5.1f ps\n', antenna(i).name, simparam(i).wn_sat);
    end
end

% close file
fclose(fid);