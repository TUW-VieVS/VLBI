% #########################################################################
% #     get_trf_data.m
% #########################################################################
%
% DESCRITPION
% This function is used to get TRF data (posistion, velocity and epoch
% info.) from different sources. The data is written to the structure
% "station". 
%   TRF data resources:
%   There are three data resources. These sources are checke for data availability in the following order.
%   If no TRF data is available in the specified TRF file (in param.txt), then the superstation file is checked, etc. 
%   The positions from position.cat serves as a kind of backup data (vel. are set to zero)
%   - 1.) TRF file (file path and name is defined in param.txt)
%   - 2.) superstation.mat (from /TRF/)
%       In this case the most recent realisation of the frame is taken. 
%       Priority order:
%           a. most recent VieTRF
%           b. vievsTRF
%           c. most recent vtrf
%           d. most recent itrf
%    - 3.) Position.cat (already available in the station struct)
%
% AUTHOR 
%   Andreas Hellerschmied, 06.12.2013
%
% INPUT
%   - PARA              : Scheduling parameter structure
%   - station           : Station structure
%
% OUTPUT
%   - station           : Station structure
%   - error_code        : Error Code (0 = no erros occured)
%   - error_msg         : Error Message (empty, if no errors occured)
%
% CHANGES
%   2015-02-12, A. Hellerschmied: added the possibility to get TRF data from superstation 
%       file and from position.cat, additionally to get it from a specific
%       TRF file (path and file name defined in PARA structure, resp. in param.txt).
%   2016-04-11, A. Hellerschmied: bug-fix: Taking coordinates from position.cat failed.
%   2017-08-31, A. Hellerschmied: Only take first 10 columns of TRF file

function [station,error_code, error_msg] = get_trf_data(station, PARA)


    % Init
    error_code = 0;
    error_msg = '';
    got_trf_data_for_station = zeros(1,length(station)); % Set field to "1", if TRF data is available for this station
    
    
    
    % ##### 1.) Try to get TRF data from TRF file 
    
    filepath_trf = PARA.TRF_FILEPATH;
    filename_trf = PARA.TRF_FILENAME; 
    
    % open TRF data file
    fid_trf = fopen([filepath_trf, filename_trf], 'r'); 
    if (fid_trf == -1)
        error_msg = 'Can not open TRF data File at the given filepath.';
        error_code = 1;
        return;
    end
    
    % Read data from TRF file:
    [temp_data] = textscan(fid_trf, '%8c %f %f %f %f %f %f %f %f %f %*[^\n]', 'CommentStyle', '%');
        
    % Get TRF data for sites specified in the "station" structure:
    number_of_stat = length(station);
    number_of_trf_datasets = length(temp_data{1}(:,1));
    
    for i_trf = 1 : number_of_trf_datasets
        stat_name_trf = temp_data{1}(i_trf,:);
        for i_stat = 1 : number_of_stat
            if (strcmp(stat_name_trf,station(i_stat).name))
                
                % Get Position
                x_trf = temp_data{2}(i_trf);
                y_trf = temp_data{3}(i_trf);
                z_trf = temp_data{4}(i_trf);
                
                % Get Velocity
                vx_trf = temp_data{5}(i_trf);
                vy_trf = temp_data{6}(i_trf);
                vz_trf = temp_data{7}(i_trf);
                
                % Get epoch, start, end
                epoch_trf = temp_data{8}(i_trf);
                start_trf = temp_data{9}(i_trf);
                end_trf = temp_data{10}(i_trf);
                
                % Save data in "station" struct
                station(i_stat).trf_xyz     = [x_trf; y_trf; z_trf];
                station(i_stat).trf_vel_xyz = [vx_trf; vy_trf; vz_trf];
                station(i_stat).trf_epoch   = epoch_trf;
                station(i_stat).trf_start   = start_trf;
                station(i_stat).trf_end     = end_trf;
                station(i_stat).trf_source  = [filepath_trf, filename_trf]; 
                
                % stat_counter = stat_counter + 1;
                got_trf_data_for_station(i_stat) = 1;
                
            end % if (strcmp(stat_name_trf,station(i_stat).name))
        end % for i_stat = 1 : number_of_stat
    end % for i_trf = 1 : number_of_trf_datasets
    
    clear temp_data;
    
    % close TRF data file
    if ( (exist('fid_trf', 'var') ) && (fid_trf ~= -1) )
        fclose(fid_trf);
    end
    
    % ##### 2.) Try to get TRF data from /TRF/superstations.mat #####
    
    % Load superstation file:
    load('../TRF/superstation.mat'); % => superstations structure
    
    %Find stations without TRF data in superstation file and get TRF data:  
    for i_stat = find(~got_trf_data_for_station)
        for i_super_stat = 1 : length(superstations)
            if strcmp(station(i_stat).name, superstations(i_super_stat).name)
                
                % select TRF from superstation file:
                % Priority order:
                % 1. most recent VieTRF
                % 2. most recent vtrf
                % 3. most recent itrf
                 
                trf_state = 'vietrf'; % Initial state
                superstat_fields = fields(superstations(i_super_stat));
                
                while(1)
                    switch(trf_state)

                        case 'vietrf'
                            if sum(strncmpi(superstat_fields, 'vietrf', 6)) > 0
                                year_max = -99;
                                trf_index = 1;
                                trf_field_names = superstat_fields(strncmpi(superstat_fields, 'vietrf', 6));
                                
                                % Run through all itrf entries:
                                for i_field = 1:length(trf_field_names)
                                    temp_field = trf_field_names{i_field};

                                    % Is field not empty?
                                    if eval(['~isempty(superstations(i_super_stat).',temp_field, ')'])
                                        % find most recent version by comparing the last 2 charscters of the field name (= year):
                                        temp_year = str2num(temp_field(length(temp_field)-1:length(temp_field)));
                                        if temp_year > year_max
                                            year_max = temp_year;
                                            trf_index = i_field;
                                        end  
                                    end
                                end
                                % Found data? 
                                if year_max > -99
                                    trf_state = 'save';
                                else
                                    trf_state = 'vievstrf';
                                end
                            else 
                                trf_state = 'vievstrf';
                            end
                            
                        case 'vievstrf'
                            if sum(strncmpi(superstat_fields, 'vievstrf', 8)) > 0
                                trf_index = 1;
                                trf_field_names = superstat_fields(strncmpi(superstat_fields, 'vievstrf', 8));
                                temp_field = trf_field_names{trf_index};
                                
                                % Is field not empty?
                                if eval(['~isempty(superstations(i_super_stat).', temp_field, ')'])
                                    trf_state = 'save';
                                else
                                    trf_state = 'vtrf';
                                end
                            else 
                                trf_state = 'vtrf';
                            end

                        case 'vtrf'
                            if sum(strncmpi(superstat_fields, 'vtrf', 4)) > 0
                                year_max = -99;
                                trf_index = 1;
                                trf_field_names = superstat_fields(strncmpi(superstat_fields, 'vtrf', 4));
                                
                                % Run through all itrf entries:
                                for i_field = 1:length(trf_field_names)
                                    temp_field = trf_field_names{i_field};

                                    % Is field not empty?
                                    if eval(['~isempty(superstations(i_super_stat).',temp_field, ')'])
                                        % find most recent version by comparing the last 2 charscters of the field name (= year):
                                        temp_year = str2num(temp_field(length(temp_field)-1:length(temp_field)));
                                        if temp_year > year_max
                                            year_max = temp_year;
                                            trf_index = i_field;
                                        end  
                                    end
                                end
                                % Found data? 
                                if year_max > -99
                                    trf_state = 'save';
                                else
                                    trf_state = 'itrf';
                                end
                            else 
                                trf_state = 'itrf';
                            end

                        case 'itrf'
                            if sum(strncmpi(superstat_fields, 'itrf', 4)) > 0
                                year_max = -99;
                                trf_index = 1;
                                trf_field_names = superstat_fields(strncmpi(superstat_fields, 'itrf', 4));
                                
                                % Run through all itrf entries:
                                for i_field = 1:length(trf_field_names)
                                    temp_field = trf_field_names{i_field};

                                    % Is field not empty?
                                    if eval(['~isempty(superstations(i_super_stat).',temp_field, ')'])
                                        % find most recent version by comparing the last 2 charscters of the field name (= year):
                                        temp_year = str2num(temp_field(length(temp_field)-1:length(temp_field)));
                                        if temp_year > year_max
                                            year_max = temp_year;
                                            trf_index = i_field;
                                        end  
                                    end
                                end
                                % Found data? 
                                if year_max > -99
                                    trf_state = 'save';
                                else
                                    trf_state = 'no_match';
                                end
                            else 
                                trf_state = 'no_match';
                            end
                        case 'no_match'
                            break;

                        case 'save'
                            
                            % Get TRF data from superstation file:
                            % Get Position
                            x_trf = eval(['superstations(i_super_stat).', trf_field_names{trf_index}, '.break.x']);
                            y_trf = eval(['superstations(i_super_stat).', trf_field_names{trf_index}, '.break.y']); 
                            z_trf = eval(['superstations(i_super_stat).', trf_field_names{trf_index}, '.break.z']);

                            % Get Velocity
                            vx_trf = eval(['superstations(i_super_stat).', trf_field_names{trf_index}, '.break.vx']);
                            vy_trf = eval(['superstations(i_super_stat).', trf_field_names{trf_index}, '.break.vy']);
                            vz_trf = eval(['superstations(i_super_stat).', trf_field_names{trf_index}, '.break.vz']);

                            % Get epoch, start, end
                            epoch_trf = eval(['superstations(i_super_stat).', trf_field_names{trf_index}, '.break.epoch']);
                            start_trf = eval(['superstations(i_super_stat).', trf_field_names{trf_index}, '.break.start']);
                            end_trf = eval(['superstations(i_super_stat).', trf_field_names{trf_index}, '.break.end']);

                            % Save data in "station" struct
                            station(i_stat).trf_xyz     = [x_trf; y_trf; z_trf];
                            station(i_stat).trf_vel_xyz = [vx_trf; vy_trf; vz_trf];
                            station(i_stat).trf_epoch   = epoch_trf;
                            station(i_stat).trf_start   = start_trf;
                            station(i_stat).trf_end     = end_trf;
                            station(i_stat).trf_source  = trf_field_names{trf_index}; 

                            got_trf_data_for_station(i_stat) = 1;
                            break;
                    end % switch(trf_state)
                end % while(1)  
            end
        end
    end
    
    
    % ##### 3.) If there is no other possibility: Take position from CATALOG file #####
    for i_stat = find(~got_trf_data_for_station) 
       % Write data in "station" struct
       station(i_stat).trf_xyz     = station(i_stat).xyz;
       station(i_stat).trf_vel_xyz = [0; 0; 0];         % default value
       station(i_stat).trf_epoch   = 51544;             % default value
       station(i_stat).trf_start   = 0;                 % default value
       station(i_stat).trf_end     = 99999;             % default value
       station(i_stat).trf_source  = 'position.cat'; 
       got_trf_data_for_station(i_stat) = 1;
    end
    
    
    % ##### 4.) Check, if there is TRF data for each station #####
    if (sum(got_trf_data_for_station) ~= length(got_trf_data_for_station))
        error_msg = ['TRF data is missing in the TRF file, superstation.mat and in position.cat for: ', station(~got_trf_data_for_station).name];
        error_code = 1;
        fclose(fid_trf);
        return; 
    end
       
return
         
    

