% #########################################################################
% #     calc_met_data
% #########################################################################
%
% DESCRIPTION
%   This function calculates meteorologic data (pressure, temperature, etc..) by evaluating GPT3.
%   The values are written to the antenna and scan structures.
%
%   #### Scan structure: ####
%   If the met data was not measured during observatuions and is missing ín the observation data 
%   (i.e. error code in NGS files of missing Met.nc file in vgosdb), GPT3 values are taken as backup!
%   They are written to the scan structure for each observations.
%   As the GPT3 values do not change rapidly with time, only one value for temp., water vapour pres. (e) and pres. is calculated at the 
%   epch of the first observation at a station (antenna(i_stat).firstObsMjd).
%
%   #### Antenna structure: ####
%   The follwoing values are calculated and store in the antenna structure (time epoch = first obs of the station; antenna(i_stat).firstObsMjd):
%    - GPT3 :
%        - pressure
%        - temperature
%        - Rel. humidity => weater vapour pressure (e)
%
%    #### Backup model and model flags: ####
%    If the backup model (GPT3) was used for AT LEAST ONE OBSERVATION of a particular antenna,
%    the according flag is set:
%        -  antenna(i_stat).gpt3temp
%        -  antenna(i_stat).gpt3pres
%        -  antenna(i_stat).gpt3e
%    These flags are also set, if the according models is defined in the GUI/parameter file.
%
%
% CREATED
%   2016-06-27     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
% - xyz2ell.m
% - gpt3_5_fast.m
%
%
% INPUT
% - parameter           - VieVS parameter structure
% - antenna             - VieVS antenna structure (from LEVEL0)
% - scan                - VieVS scan structure (from LEVEL0)
%
%
% OUTPUT
% - antenna             - VieVS antenna structure (updated with met. data)
% - scan                - VieVS scan structure (updated with met. data)
%
% CHANGES:
% 2016-07-06, A. Hellerschmied: - Minor bug fixed (check of e in the scan structure)
%                                 - Flags in antenna struct (gpt2pres, gpt2temp, gpt2e, gpttemp, gptpres) are no set:
%                                    => They can be used e.g. to check after calculating the met. data, if GPT or GPT2 values were used as backup! 
% 2017-01-23, D. Landskron: - all GPT2 changed to GPT3, GPT removed
% 2017-04-04, D. Landskron: bug corrected with scan.stat.e GPT3 backup handling
%

function [antenna, scan] = calc_met_data(parameter, antenna, scan)

% Init:
number_of_stations  = length(antenna);
number_of_scans     = length(scan);

% Preallocation:
gpt3_p      = zeros(number_of_stations, 1);
gpt3_T      = zeros(number_of_stations, 1);
gpt3_dT     = zeros(number_of_stations, 1);
gpt3_e      = zeros(number_of_stations, 1);
gpt3_ah     = zeros(number_of_stations, 1);
gpt3_aw     = zeros(number_of_stations, 1);
gpt3_undu   = zeros(number_of_stations, 1);
gpt3_Tm     = zeros(number_of_stations, 1);
gpt3_lambda = zeros(number_of_stations, 1);
gpt3_Gn_h   = zeros(number_of_stations, 1);
gpt3_Ge_h   = zeros(number_of_stations, 1);
gpt3_Gn_w   = zeros(number_of_stations, 1);
gpt3_Ge_w   = zeros(number_of_stations, 1);

% Options
error_code_invalid_met_data = -999; % Error corde for missing met. data in NGS file (numerical), set in read_ngs.m


% ##### ANTENNA #####

% Preallocation:

% read the cell_grid of GPT3
cell_grid_GPT3 = gpt3_5_fast_readGrid;

% #### Loop over all antennas ####
for i_stat = 1 : number_of_stations
    
    % Get time epoch for further calculations:
    mjd = antenna(i_stat).firstObsMjd;
    
    % Ell. coordinates:
    [lat, lon, h] = xyz2ell([antenna(i_stat).x, antenna(i_stat).y, antenna(i_stat).z]);
    
    
    % #### Calc. GPT3 data (call gpt3_5_fast.m) ####
    [gpt3_p(i_stat) , gpt3_T(i_stat) , gpt3_dT(i_stat) , gpt3_Tm(i_stat) , gpt3_e(i_stat) , gpt3_ah(i_stat) , gpt3_aw(i_stat), gpt3_lambda(i_stat) , gpt3_undu(i_stat) , gpt3_Gn_h(i_stat) , gpt3_Ge_h(i_stat) , gpt3_Gn_w(i_stat) , gpt3_Ge_w(i_stat) ] = gpt3_5_fast (mjd,lat,lon,h,0,cell_grid_GPT3);
    
    % Save GPT3 data to antenna struct.:
    antenna(i_stat).gpt3.p      = gpt3_p(i_stat);
    antenna(i_stat).gpt3.T      = gpt3_T(i_stat);
    antenna(i_stat).gpt3.dT     = gpt3_dT(i_stat);
    antenna(i_stat).gpt3.Tm     = gpt3_Tm(i_stat);
    antenna(i_stat).gpt3.e      = gpt3_e(i_stat);
    antenna(i_stat).gpt3.ah     = gpt3_ah(i_stat);
    antenna(i_stat).gpt3.aw     = gpt3_aw(i_stat);
    antenna(i_stat).gpt3.lambda = gpt3_lambda(i_stat);
    antenna(i_stat).gpt3.undu   = gpt3_undu(i_stat);
    antenna(i_stat).gpt3.Gn_h   = gpt3_Gn_h(i_stat);
    antenna(i_stat).gpt3.Ge_h   = gpt3_Ge_h(i_stat);
    antenna(i_stat).gpt3.Gn_w   = gpt3_Gn_w(i_stat);
    antenna(i_stat).gpt3.Ge_w   = gpt3_Ge_w(i_stat);
    

    % #### Set flags ####
    % => If GPT3 is defined in the GUI/parameters structure
    % Temp.:
%     if strcmp(parameter.vie_init.tp,'gpt3')
%         antenna(i_stat).gpt3temp = 1;
%     end
%     % Pres.:
%     if strcmp(parameter.vie_init.pt,'gpt3')
%         antenna(i_stat).gpt3pres = 1;
%     end
%     % Hum.:
%     if strcmp(parameter.vie_init.pt,'gpt3')
%         antenna(i_stat).gpt3e = 1;
%     end
    
end


% ##### SCAN #####

% #### Loop over all scans ####
for i_scan = 1 : number_of_scans
   
    % #### Loop over all stations ####
    % Each scan has one field for all stations in the antenna structure
    %   => Fill them with met. data, even if a particular station does not contribute to a scan
    max_number_of_stations_in_scan = length(scan(i_scan).stat);
    
    for i_stat = 1 : max_number_of_stations_in_scan
        
        % Temp.:
        if ~isempty(scan(i_scan).stat(i_stat).temp)
            if ~((scan(i_scan).stat(i_stat).temp ~= error_code_invalid_met_data)  &&  strcmp(parameter.vie_init.tp,'in situ'))
                switch(parameter.vie_init.tp)
                    case {'in situ'}
                        scan(i_scan).stat(i_stat).temp  = gpt3_T(i_stat);
                        antenna(i_stat).gpt3temp        = 1;
                    case {'gpt3'}
                        scan(i_scan).stat(i_stat).temp  = gpt3_T(i_stat);
                end
            end
        end
        
        % Pres.:
        if ~isempty(scan(i_scan).stat(i_stat).pres)
            if ~((scan(i_scan).stat(i_stat).pres  ~= error_code_invalid_met_data)  &&  strcmp(parameter.vie_init.zhd,'in situ'))
                switch(parameter.vie_init.zhd)
                    case {'in situ'}
                        scan(i_scan).stat(i_stat).pres  = gpt3_p(i_stat);
                        antenna(i_stat).gpt3pres        = 1;
                    case {'gpt3'}
                        scan(i_scan).stat(i_stat).pres  = gpt3_p(i_stat);
                end
            end
        end
        
        % e:
        if ~isempty(scan(i_scan).stat(i_stat).e)
            if ~((scan(i_scan).stat(i_stat).e  ~= error_code_invalid_met_data)  &&  strcmp(parameter.vie_init.zwd,'in situ'))
                switch(parameter.vie_init.zwd)
                    case {'in situ'}
                        scan(i_scan).stat(i_stat).e = gpt3_e(i_stat);
                        antenna(i_stat).gpt3e       = 1;
                    case {'gpt3'}
                        scan(i_scan).stat(i_stat).e = gpt3_e(i_stat);
                        antenna(i_stat).gpt3e       = 1;
                    case {'no'}
                        if scan(i_scan).stat(i_stat).e == error_code_invalid_met_data
                            scan(i_scan).stat(i_stat).e = gpt3_e(i_stat);
                            antenna(i_stat).gpt3e       = 1;
                        end
                end
            end
        end
        
    end % for i_stat = 1 : number_of_stations
            
end % for i_scan = 1 : number_of_scans
            

    
