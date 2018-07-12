% #########################################################################
% #     read_vgosdb_input_settings
% #########################################################################
%
% DESCRIPTION
% This function reads the vgosDB wrapper files and strores the information
% in a MATLAB structure.
% 
% Attention: Matlab does not allow to use '-' in field names. Therefore, all
%            according '-' (e.g. in stationnames) are replaced by '_'.
%            E.g.: "BR-VLBA" => "BR_VLBA"
%
% CREATED  
%   2018-06-12     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%
%
% INPUT:
%   - path_nc
%   - session_name
%	- institute (institute name)
%   - wrapper_k  (wrapper tag, "all" or "ngs")
%   - wrapper_v  (wrapper version number)
% OUTPUT:
%	- wrapper_data
%
% CHANGES:
%

function [ wrapper_data ] = read_vgosdb_wrapper(path_nc, session_name, institute, wrapper_k, wrapper_version)

% ##### Options #####

% Wrapper blocks which are considered (all others are neglected...):
wraper_blocks = {'Session', 'Station', 'Scan', 'Observation'};
wrapper_data = struct('Session', [], 'Station', [], 'Scan', [], 'Observation', [], 'wrapper_filename', []);


% ##### Get wrapper filename #####

% Get wrapper filenames in directory:
tmp = dir(path_nc);
tmp = {tmp.name};
wrapper_filenames = {};
i_wrapper = 0;
version = 0;
for i = 1 : length(tmp)
    name = tmp{i};
    if length(name) > 5
        if strcmp(name(end-3:end), '.wrp')
            i_wrapper = i_wrapper + 1;
            wrapper_filenames{i_wrapper} = name;
            
            % version number
            str_match = ['_i', institute, '_k', wrapper_k, '.wrp'];
            if strcmp(name(end-(length(str_match)-1):end), str_match)
                id = strfind(name, '_V');
                if str2num(name(id+2:id+4)) > version
                    version = str2num(name(id+2:id+4));
                end
            end
        end
    end
end



% Get version number
switch(wrapper_version)
    case 'highest_version'
        version_str = sprintf('%03d', version);
        
    otherwise
        version_str = sprintf('%03d', str2num(wrapper_version));
end

% disp(version_str)
nc_filename = [session_name, '_V', version_str, '_i', institute, '_k', wrapper_k, '.wrp'];
wrapper_data.wrapper_filename = nc_filename;
    

% ##### Open file #####
fid = fopen([path_nc, nc_filename],'r');

if fid<0
    error('%s does not exist\n', [path_nc, nc_filename])
else
    fprintf('Reading wrapper: %s\n', [path_nc, nc_filename])
    
    % ##### Loop over all lines #####
    
    % Loop init:
    flag_in_block = 0;
    flag_valid_block = 0;
    current_block_str   = '';
    current_default_dir = '';
    
    
    while ~feof(fid)
	
        % #### Read one line ####
        str = fgetl(fid);
        % disp(str);
        
        if isempty(str)
            continue;
        end 
        if strcmp(str(1), '!') % comment
            continue;
        end
        temp_str = textscan(str, '%s'); 
        
        switch(flag_in_block)
            
            case 0
                % Begin of block?
                if strcmp(temp_str{1}{1}, 'Begin')
                    flag_in_block = 1;
                    current_block_str = temp_str{1}{2};
                    if sum(strcmp(temp_str{1}{2}, wraper_blocks)) > 0 % Begin of valid wrapper block found
                        flag_valid_block = 1;
                    else
                        flag_valid_block = 0;
                    end
                else
                    continue;
                end
                
            case 1 % in block
                switch(flag_valid_block)
                    
                    case 0 % NOT IN valid block => wait for end of this block!
                        if strcmp(temp_str{1}{1}, 'End') && strcmp(temp_str{1}{2}, current_block_str)  % End of block?
                            flag_in_block = 0;
                        end
                       
                    case 1 % IN valid block
                        if strcmp(temp_str{1}{1}, 'End') && strcmp(temp_str{1}{2}, current_block_str)  % End of block?
                            flag_in_block = 0;
                            current_default_dir = '';
                            continue;
                        end
                        if strcmp(temp_str{1}{1}, 'Default_Dir')
                            current_default_dir = temp_str{1}{2};
                            % Exchange '-' with '_', because '-' are not valid for field names of MATLAB structs:
                            current_default_dir(strfind(current_default_dir, '-')) = '_';
                            continue;
                        end
                        % nc file?
                        if size(temp_str{1}, 1) == 1
                            tmp = temp_str{1}{1};
                            if strcmp(tmp(end-2:end), '.nc') && ~isempty(current_default_dir) % nc file! => make entry!
                                if ~isfield(wrapper_data.(current_block_str), current_default_dir)
                                    wrapper_data.(current_block_str).(current_default_dir).files = temp_str{1}{1};
                                else
                                    wrapper_data.(current_block_str).(current_default_dir).files = [wrapper_data.(current_block_str).(current_default_dir).files, temp_str{1}];
                                end
                            end
                        end
                         
                        
                end % switch(flag_valid_block)
            
        end % switch(flag_in_valid_block)
            
        
    end % while ~feof(fid)
    
    
    
    
    % Close file
    fclose(fid);

    
end




