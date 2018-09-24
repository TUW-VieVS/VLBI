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
% 2018-08-28, D. Landskron: adapted also for different wrapper file name
%
% #########################################################################



function [ wrapper_data ] = read_vgosdb_wrapper(path_nc, session_name, institute, wrapper_k, wrapper_version)

% ##### Options #####

% Wrapper blocks which are considered (all others are neglected...):
wraper_blocks = {'Session', 'Station', 'Scan', 'Observation'};
wrapper_data = struct('Session', [], 'Station', [], 'Scan', [], 'Observation', [], 'wrapper_filename', []);
min_wrapper_version = 4;


% #### Init #####
num_institutions    = length(institute);


% ##### Get wrapper filename #####
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
           
        end
    end
end

num_of_length_wrappers = length(wrapper_filenames);

% Get i, k, V flags of wrapper names
i_list = cell(num_of_length_wrappers, 1);
V_list = nan(num_of_length_wrappers, 1);
k_list = cell(num_of_length_wrappers, 1);

for i = 1 : num_of_length_wrappers
    % _i
    str_id_start = strfind(wrapper_filenames{i}, '_i');
    if  ~isempty(str_id_start)
        substr = wrapper_filenames{i}(str_id_start+2:end-4);
        tmp_id = strfind(substr, '_');
        if isempty(tmp_id)
            i_list{i} = substr;
        else
            str_id_end = min(tmp_id);
            i_list{i} = substr(1:str_id_end-1);
        end
    end
    
    % _k
    str_id_start = strfind(wrapper_filenames{i}, '_k');
    if  ~isempty(str_id_start)
        substr = wrapper_filenames{i}(str_id_start+2:end-4);
        tmp_id = strfind(substr, '_');
        if isempty(tmp_id)
            k_list{i} = substr;
        else
            str_id_end = min(tmp_id);
            k_list{i} = substr(1:str_id_end-1);
        end
    end
    
    % _V
    str_id_start = strfind(wrapper_filenames{i}, '_V');
    if  ~isempty(str_id_start)
        substr = wrapper_filenames{i}(str_id_start+2:end-4);
        tmp_id = strfind(substr, '_');
        if isempty(tmp_id)
            V_list(i) = str2num(substr);
        else
            str_id_end = min(tmp_id);
            V_list(i) = str2num(substr(1:str_id_end-1));
        end
    end

end

% flag arraycontent:
%  col 1: version
%  col 2: k
%  col 3 to n: institution(s)
flag_array = zeros(num_of_length_wrappers, (2 + num_institutions));

nc_filename = '';


% ##### Select a wrapper file by the following conditions: #####
% - Version number:
%   - if there is no specfic version number defined, the highes number is taken of the defined institution (first in the list)
%     - the min wrapper version is defined in the variable "min_wrapper_version"
%   - if a version number is defined, the wrapper with this version and the institution with highest priority is taken

% _V
if isempty(wrapper_version)
    % any k is valid!
    flag_array(V_list >= min_wrapper_version, 1) = V_list(V_list >= min_wrapper_version);
else
    if strcmp(wrapper_version, 'highest_version')
        flag_array(V_list >= min_wrapper_version, 1) = V_list(V_list >= min_wrapper_version);
    else % only take specified version:
       flag_array(:, 1) = V_list == wrapper_version;
    end
end

% _k
if isempty(wrapper_k)
    % any k is valid!
    flag_array(:, 2) = true;
else
    flag_array(:, 2) = strcmp(k_list, 'all');
end

%_i
for i_int = 1 : num_institutions
    i_col = i_int+2;
    flag_array(:, i_col) = strcmp(i_list, institute{i_int});
end


% Select the best fitting wrapper file:

% exclude filenames based on _k
wrapper_filenames_tmp = wrapper_filenames(logical(flag_array(:, 2)));
flag_array_tmp =  flag_array(logical(flag_array(:, 2)), :);
if isempty(wrapper_filenames_tmp)
    error(' - No valid wrapper file found (no wrapper with _k = %s)',wrapper_k);
end

% Loop over institutions and select the highest version available:
for i_inst = 1 : num_institutions
   i_col = i_inst + 2;

   if sum(flag_array_tmp(:, i_col)) > 0 % Check, if there are valid wrappers with this inst.
       % Check the version numbers for this inst.:
       if sum(flag_array_tmp(:,1) & flag_array_tmp(:, i_col)) == 0
           % error('No valid wrapper file found!)');
           fprintf(' - No valid wrapper file for institution: %s\n', institute{i_inst})
           continue
       end
       flag_array_tmp(~(flag_array_tmp(:,1) & flag_array_tmp(:, i_col)), 1) = 0;
       [max_version, max_version_id]=max(flag_array_tmp(:, 1));
       
       nc_filename = wrapper_filenames_tmp{max_version_id};
   else
       fprintf(' - No valid wrapper file for institution: %s\n', institute{i_inst})
       
   end
   
end

% Check, if wrapper is available!
if isempty(nc_filename)
    error('No valid wrapper file found!');
end

% Load wrapper file:
fid = fopen([path_nc, nc_filename],'r');
if fid < 0
    error('Unable to open wrapper file: %s', [path_nc, nc_filename]);
end


wrapper_data.wrapper_filename = nc_filename;

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




