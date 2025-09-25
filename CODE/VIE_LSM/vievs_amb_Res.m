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
%
% ************************************************************************

function [amb2] = vievs_amb_Res(parameter, antenna, sources, scan, res)

tic;

amb2 = 1;
fprintf('\n');
fprintf('Ambiguities Residual Correction: On \n');

% ambSpace in scan-File
ambspace = (scan(1).obs(1).ambspace).* 10^9; % Ambiguity Spacing [ns]
ambspaceCM = abs(ambspace*30); % Ambiguity Spacing positive [cm]
% ambspaceCM1 = ambspaceCM*0.7;
% ambspaceCM2 = ambspaceCM*1.3;
% Border values because ambiguity value does not match with exact value
ambspaceNS1 = ambspace*0.8; % lower border
ambspaceNS2 = ambspace*1.2; % upper border

resFirst = res.firstVal; % First Solution Residuals [cm]
AmbValue = zeros(size(resFirst));

for k = 1:length(resFirst)
    if abs(resFirst(k)) > ambspaceCM/2
        diffs = resFirst(k)-resFirst;
        for jj = 1:length(resFirst) % try to find remaining ambiguities in the residuals
            diffs_ns = diffs(jj)/30;
            if sign(ambspaceNS1) == sign(diffs_ns) % positive ambiguities
                if diffs_ns > ambspaceNS1 && diffs_ns < ambspaceNS2
                    AmbValue(k) = ambspaceCM;
                    break;
                elseif diffs_ns > 2*ambspaceNS1 && diffs_ns < 2*ambspaceNS2
                    AmbValue(k) = 2*ambspaceCM;
                    break;
                end
            elseif sign(-ambspaceNS1) == sign(diffs_ns) % negative ambiguities
                if diffs_ns > -ambspaceNS2 && diffs_ns < -ambspaceNS1
                    AmbValue(k) = -ambspaceCM;
                    break;
                elseif diffs_ns > -2*ambspaceNS2 && diffs_ns < -2*ambspaceNS1
                    AmbValue(k) = -2*ambspaceCM;
                    break;
                end
            end
        end
    end

    % if sign(ambspaceCM1) == sign(resFirst(k))
    %     if resFirst(k) > ambspaceCM1 && resFirst(k) < ambspaceCM2
    %         AmbValue(k) = ambspaceCM;
    %         % res.firstVal(k) = resFirst(k) + ambspaceCM;
    %     end
    % elseif sign(-ambspaceCM1) == sign(resFirst(k))
    %     if resFirst(k) > -ambspaceCM2 && resFirst(k) < -ambspaceCM1
    %         AmbValue(k) = -ambspaceCM;
    %         % res.firstVal(k) = resFirst(k) - ambspaceCM;
    %     end
    % end
end


% AmbValue
% save('AmbValue.mat','AmbValue')



%% Save Ambiguities to TXT File

% Ambiguity file
parameter.amb.amb_file_dir = 'PU';
parameter.amb.flag_change_amb = true; %true

checkPath = ['../DATA/AMB/', parameter.amb.amb_file_dir, '/', parameter.filepath(end-4:end-1), '/'];
if ~exist(checkPath,'dir')
    mkdir(checkPath);
end

output_file_path = ['../DATA/AMB/', parameter.amb.amb_file_dir, '/', parameter.filepath(end-4:end-1), '/', [parameter.session_name '_' parameter.vie_init.vgosDb_observation_parameter(end) ], '.AMB'];

exist(output_file_path, 'file');

if parameter.amb.flag_change_amb
    if exist(output_file_path, 'file')
        % [parameter.amb.obs2change] = readAMB(amb_filename_path);
        % fprintf('%d baselines with ambiguities will be changed\n',size(parameter.amb.obs2change,2)); 
        % scan = changeAMB(scan, antenna, sources, parameter);
    else
        fprintf('Ambiguity list not available: %s\n', output_file_path);
    end
else
    fprintf('Ambiguities will not be changed\n');
end
% fprintf('\n')

file_id = fopen(output_file_path,'a+'); % CHECK AND MAYBE CHANGE!!!

ambNumberCount = 0;

% Iterate through all rows and columns of the adjustmentsTable
for k = 1:length(resFirst)

    % station 1
    station1_num = res.baselineOfObs(k,1);
    station1_name = strtrim(antenna(station1_num).name); 
    % station 2
    station2_num = res.baselineOfObs(k,2);
    station2_name = strtrim(antenna(station2_num).name);
    % mjd
    mjd = res.mjd(k);
    % iso
    iso_idx = res.source(k); % ISO index for source name
    % source name with iso nummer
    source_name = strtrim(sources.q(iso_idx).name); % no extra spaces
    % amb value
    amb_value1 = AmbValue(k)/30 ; % [ns]
    amb_value = round(amb_value1);

    if amb_value == 0 % Skip if ambiguity is zero
        continue;
    end

    formatted_mjd = sprintf('%.12f', mjd); % Format MJD with 12 decimal places

    if amb_value > 0 % Remove leading zeros for ambiguity and format
        formatted_corr = sprintf('%d', round(amb_value)); % Positive value
    else
        formatted_corr = sprintf('-%d', abs(round(amb_value)));
    end

    % Format the output line in the required fixed-width format
    line_data = {sprintf('%-8s', station1_name), sprintf('%-8s', station2_name), ...
                sprintf('%-18s', formatted_mjd), sprintf('%-8s', source_name), ...
                sprintf('%-5s', formatted_corr)};

    ambNumberCount = ambNumberCount + 1;   
    fprintf('\t \t %s\n', strjoin(line_data, ' '));
            
    % Write the line to the file, using spaces for separation
    fprintf(file_id, '%s\n', strjoin(line_data, ' ')); 
    
end

fclose(file_id);
fprintf('\t A total of %d ambiguities were found. \n', ambNumberCount);
fprintf('\t Ambiguities written to %s\n', output_file_path);
fprintf('\t After successful calculation, please deselect the compute checkbox, select the apply checkbox and rerun vie_lsm.\n');

endtime = toc;

fprintf('Ambiguities Residual Correction: finished after %.1f seconds \n', endtime')

% Count the number of NaN values and 0
% numNaNs = sum(isnan(diffs));
% numZeros = sum(diffs == 0);
% disp(['Number of NaN values: ', num2str(numNaNs)]);
% disp(['Number of zero values: ', num2str(numZeros)]);

end