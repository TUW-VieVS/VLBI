% ************************************************************************
%   Description:
%   This function is used after the automatic ambiguity resolution.
%   Remaining ambiguities are calculated based on the first solution.
%
%   Input:	
%      Variables from parameter, antenna, sources, scan, res
%
%   Output:
%      Ambiguity values for group delays saved into .AMB file 
%
%   Coded for VieVS: 
%       25.09.2025 by Peter Urban
%
%   Update:
%       15.10.2025 by Peter Urban: several smaller improvements
%       23.04.2026 by Peter Urban: major update of the approach for
%       handling different ambiguity spacings per baselines, small bug
%       fixes
%
% ************************************************************************

function [amb2] = vievs_amb_Res(parameter, antenna, sources, scan, res)

tic;

amb2 = 1;
fprintf('\n');
fprintf('Ambiguities Residual Correction: On \n');

idx_X = 1;
ambspaceVec = [];

% data read
for i = 1:length(scan) 
    for j = 1:length(scan(i).obs) 
        ambspaceVec(idx_X) = (scan(i).obs(j).ambspace).* 10^9;
        idx_X = idx_X + 1;
    end
end

% find Unique ambspace and remove 0 and NaN and sort descending
ambspace_unique = unique(ambspaceVec);
ambspace_unique(ambspace_unique == 0 | isnan(ambspace_unique)) = []; 
ambspace_unique = sort(ambspace_unique, 'descend'); 


resFirst = res.firstVal; % First Solution Residuals [cm]
AmbValue = zeros(size(resFirst)); 

% main loop for ambspace_unique
for uIdx = 1:length(ambspace_unique)
    
    % set ambspace
    curr_ambspace = ambspace_unique(uIdx);
    
    fprintf('Check ambspace: %.4f ns ... ', curr_ambspace);
    
    ambspaceCM = abs(curr_ambspace * 30);   % Ambiguity Spacing [cm]
    
    % borders
    base_lower = curr_ambspace * 0.80; % for X-Band: 0.9
    base_upper = curr_ambspace * 1.20; % for X-Band: 1.1
    
    % Reset AmbValue 
    AmbValue = zeros(size(resFirst));
    
    % Calc
    for k = 1:length(resFirst)
        % check if it is too big and needs maybe correction
        if abs(resFirst(k)) > ambspaceCM/2
            
            diffs = resFirst(k) - resFirst;
            
            for jj = 1:length(resFirst) % check with other residuals
                diffs_ns = diffs(jj) / 30;
                
                % check for +-10x ambspace
                found_match = false;
                for m = 1:10
                    lower_lim = m * base_lower;
                    upper_lim = m * base_upper;
                    
                    if sign(curr_ambspace) == sign(diffs_ns) % Positive Ambiguities
                        if diffs_ns > lower_lim && diffs_ns < upper_lim
                            AmbValue(k) = m * ambspaceCM;
                            found_match = true;
                        end
                    elseif sign(-curr_ambspace) == sign(diffs_ns) % Negative Ambiguities
                        if diffs_ns > -upper_lim && diffs_ns < -lower_lim
                            AmbValue(k) = -m * ambspaceCM;
                            found_match = true;
                        end
                    end
                    
                    if found_match
                        break;
                    end
                end
                
                if found_match
                    break; 
                end
            end
        end
    end
    
    % Check if solution got found
    current_sum = sum(abs(AmbValue));
    
    if current_sum > 0
        fprintf('Success! Ambiguities found with ambspacing %.4f ns.\n', curr_ambspace);
        break; 
    else
        fprintf('No corrections found. Try next ambspacing. \n');
    end
    
end

if sum(abs(AmbValue)) == 0
    fprintf('Warning: No ambiguities with ambspace_unique found.\n');
end



%% Save Ambiguities to TXT File

% Ambiguity file
parameter.amb.flag_change_amb = true; %true

if strcmp(parameter.vie_init.AmbDir, parameter.vie_init.level0OutDir)
    AmbDir = parameter.vie_init.level0OutDir;
elseif isempty(parameter.vie_init.AmbDir)
    AmbDir = parameter.vie_init.level0OutDir;
else
    AmbDir = parameter.vie_init.AmbDir;
end

checkPath = ['../DATA/AMB/', AmbDir, '/'];
if ~exist(checkPath,'dir')
    mkdir(checkPath);
end

output_file_path = ['../DATA/AMB/', AmbDir, '/', [parameter.session_name '_' parameter.vie_init.vgosDb_observation_parameter(end) ], '.AMB'];

exist(output_file_path, 'file');

if parameter.amb.flag_change_amb
    if exist(output_file_path, 'file')
        fid = fopen(output_file_path,'r');
        a = 1;
        while ~feof(fid)
            str = fgetl(fid);
            obsamb(a).sta1 = str(1:8);
            obsamb(a).sta2 = str(10:17);
            obsamb(a).mjd = str2num(str(19:36));
            obsamb(a).sou = str(38:45);
            obsamb(a).amb = str2num(str(47:end));
            a=a+1;
        end
        fclose(fid);
        OrigObsLen = length([obsamb.mjd]');
    else
        fprintf('Ambiguity list not available: %s\n', output_file_path);
    end
else
    fprintf('Ambiguities will not be changed\n');
end

file_id = fopen(output_file_path,'w+'); % CHECK AND MAYBE CHANGE!!!
ambNumberCount = 0;

% Iterate through all rows and columns of the adjustmentsTable
for k = 1:OrigObsLen

    station1_name = strtrim(obsamb(k).sta1); 
    station2_name = strtrim(obsamb(k).sta2);
    mjd = obsamb(k).mjd;
    source_name = strtrim(obsamb(k).sou); % no extra spaces
    ambvalueFile = obsamb(k).amb;
    
    
    countttt = 0;
    amb_value2 = 0;
    
    for kk = 1:length(resFirst)
        % station 1
        station1_numR = res.baselineOfObs(kk,1);
        station1_nameR = strtrim(antenna(station1_numR).name); 
        % station 2
        station2_numR = res.baselineOfObs(kk,2);
        station2_nameR = strtrim(antenna(station2_numR).name);
        % mjd
        mjdR = res.mjd(kk);
        % iso
        iso_idxR = res.source(kk); % ISO index for source name
        % source name with iso nummer
        source_nameR = strtrim(sources.q(iso_idxR).name); % no extra spaces
        if mjd == mjdR
            if strcmp(station1_name, station1_nameR) && strcmp(station2_name, station2_nameR) && strcmp(source_name, source_nameR)
                amb_value2 = (AmbValue(kk)/30); % [ns]
                countttt = countttt +1;   
                if amb_value2 ~= 0
                     formatted_mjd1 = sprintf('%.12f', mjd); % Format MJD with 12 decimal places
                    if amb_value2 > 0
                        formatted_corr1 = sprintf('%d', round(amb_value2));
                    elseif amb_value2 == 0
                        formatted_corr1 = sprintf('%d', round(amb_value2));
                    else
                        formatted_corr1 = sprintf('-%d', abs(round(amb_value2)));
                    end
                        % Format the output line in the required fixed-width format
                    line_data1 = {sprintf('%-8s', station1_name), sprintf('%-8s', station2_name), ...
                    sprintf('%-18s', formatted_mjd1), sprintf('%-8s', source_name), ...
                    sprintf('%-5s', formatted_corr1)};
                    fprintf('\t \t %s\n', strjoin(line_data1, ' '));
                end
            end
        end
    end
    
    % amb value
    amb_value1 = amb_value2 + ambvalueFile; % [ns]
    amb_value = round(amb_value1);

    formatted_mjd = sprintf('%.12f', mjd); % Format MJD with 12 decimal places

    if amb_value > 0
            formatted_corr = sprintf('%d', round(amb_value));
    elseif amb_value == 0
            formatted_corr = sprintf('%d', round(amb_value));
    else
            formatted_corr = sprintf('-%d', abs(round(amb_value)));
    end

    % Format the output line in the required fixed-width format
    line_data = {sprintf('%-8s', station1_name), sprintf('%-8s', station2_name), ...
                sprintf('%-18s', formatted_mjd), sprintf('%-8s', source_name), ...
                sprintf('%-5s', formatted_corr)};

    ambNumberCount = ambNumberCount + 1;   

    % Write the line to the file, using spaces for separation
    fprintf(file_id, '%s\n', strjoin(line_data, ' ')); 

end


end