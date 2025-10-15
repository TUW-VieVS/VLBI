% ************************************************************************
%   Description:
%   This function calculates automatically the ambiguity resolution based 
%   on multiband-singleband delay differences of the group delays with 
%   triangle delay closures.
%
%   Input:	
%      Variables from out_struct, wrapper_data, parameter, antenna,
%      sources, scan
%
%   Output:
%      Ambiguity values for group delays saved into .AMB file 
%
%   Coded for VieVS: 
%       25.09.2025 by Peter Urban
%
%   Update:
%       15.10.2025 by Peter Urban: adjustment of the threshold maximum
%       closure limit and several smaller improvements
%
% ************************************************************************

function [amb1] = vievs_amb(out_struct,wrapper_data, parameter, antenna, sources, scan)


%% Load and prepare variables

tic;

format longg

fprintf('\t Ambig corr: directly calculated in VieVS \n')

amb1 = 1;

ambspace = (out_struct.Observables.AmbigSize_bX.AmbigSize.val).* 10^9; % Ambiguity spacing [ns]
% ambspace = min(ambspace);
ambspace = mode(ambspace);
threshold = ambspace/2;  % Threshold for closure value
% threshold = 10;  % Threshold for closure value

medianSpecial = 0; % check if a special case for median calculation exist

totalObservations=out_struct.head.NumObs.val; % Number of observations

% Multiband and Singleband Group Delay for X-Band
SBD_X = double(out_struct.Observables.SBDelay_bX.SBDelay.val);
MBD_X = double(out_struct.Observables.GroupDelay_bX.GroupDelay.val);

diffX = (MBD_X - SBD_X).* 10^9; % mb-sb differences [ns]

idx_X = 1;
mb_sb_X(totalObservations) = struct('scan', [], 'obs', [], 'diff', [], 'i1', [], 'i2', [], 'mjd', [], 'isoV', []);
mnn = length(antenna); % Total number of stations

% Read all the data and save it for calculations
for i = 1:length(scan) 
    for j = 1:length(scan(i).obs) 
        mb_sb_X(idx_X).scan = i;
        mb_sb_X(idx_X).obs = j;
        mb_sb_X(idx_X).diff = diffX(idx_X); % [ns]
        mb_sb_X(idx_X).i1 = scan(i).obs(j).i1;
        mb_sb_X(idx_X).i2 = scan(i).obs(j).i2;
        mb_sb_X(idx_X).mjd = scan(i).mjd;
        mb_sb_X(idx_X).isoV = scan(i).iso;
        idx_X = idx_X + 1;
        % if i == 4 || j == 4 % Exclude Station 3 for one particular Session
        %     mb_sb_X(idx_X).diff = NaN; 
        % end
    end
end



%% Ambiguity Calculation with Iterative Baseline Correction

allBaselines = NaN(mnn, mnn);  % Table for all baseline medians
combinations = nchoosek(1:mnn, 3); % Generate all combinations of 3 stations 

% comb_orig = combinations;
% comb_rand = combinations(randperm(size(combinations,1)), :);
% start_row = 86;
% comb_shifted = [combinations(start_row:end, :);
%     combinations(1:start_row-1, :)];
% save('comb_shifted.mat');
% idxx = find(combinations(:,1) == 4, 1,"first");
% if ~isempty(idxx)
%     comb_shifted = [combinations(idxx:end, :);
%         combinations(1:idxx-1, :)];
% end

i1Indices = [mb_sb_X.i1]';
i2Indices = [mb_sb_X.i2]';
diffs = [mb_sb_X.diff]';

for i = 1:mnn
    for j = 1:mnn
        if i == j, continue; end  % Skip diagonal entries

        % Find all data for baseline i-j
        indices = find(i1Indices == i & i2Indices == j);
        if isempty(indices)
            continue;  % No data for this baseline
        end

        diffsForPair = diffs(indices);
        adjustedMedian = calculateAdjustedMedian(diffsForPair, ambspace);
        allBaselines(i, j) = adjustedMedian;
    end
end

% Table for the storage of ambiguity values
adjustmentsTable = array2table(allBaselines, 'VariableNames', arrayfun(@(x) sprintf('S%d', x), 1:mnn, 'UniformOutput', false));

% Replace NaN values in adjustmentsTable with NaN and all non-NaN values with 0
for row = 1:mnn
    for col = 1:mnn
        if isnan(allBaselines(row, col))
            if row < col && ~isnan(allBaselines(col, row))
                % If in the upper triangle (row < col) NaN, but there is a value in lower triangle (col, row) 
                adjustmentsTable{row, col} = 0;
                allBaselines(row, col) = -allBaselines(col, row);
            else
                adjustmentsTable{row, col} = NaN; % keep NaN-value
            end
        else
            adjustmentsTable{row, col} = 0; % exchange not-NaN-values through 0
        end
    end
end

iteration_changes = true; % Flag to track if any adjustment was made in an iteration
iteration_count = 0;  % Counter for number of iterations
total_iteration_count = 0;

while iteration_changes
    iteration_count = iteration_count + 1;
    total_iteration_count = total_iteration_count + 1;
    iteration_changes = false;  % Reset flag at the beginning of each iteration

    for idx = 1:size(combinations, 1)
        i = combinations(idx, 1);  % Station 1 (i1)
        j = combinations(idx, 2);  % Station 2 (i2)
        k = combinations(idx, 3);  % Station 3 (i3)

        % Compute the Closure value only if there are no NaN values
        median_diff21 = allBaselines(i, j);  % i1 - i2
        median_diff32 = allBaselines(j, k);  % i2 - i3
        median_diff13 = allBaselines(i, k);  % i1 - i3

        % Check for NaN values and skip calculation if any are present
        if isnan(median_diff21) || isnan(median_diff32) || isnan(median_diff13)
            % fprintf('   Triangle (i=%d, j=%d, k=%d):\n', i, j, k);
            % fprintf('     Baseline values: %.2f, %.2f, %.2f\n', median_diff21, median_diff32, median_diff13);
            % fprintf('     Skipping calculation due to NaN value(s).\n');
            continue;  % Skip this combination if any value is NaN
        end

        % Centroid calculation
        P1 = [antenna(i).x, antenna(i).y, antenna(i).z]';
        P2 = [antenna(j).x, antenna(j).y, antenna(j).z]';
        P3 = [antenna(k).x, antenna(k).y, antenna(k).z]';
        sv = (P1 + P2 + P3) / 3;  % Centroid

        % Basis vectors from the stations
        b13 = P3 - P1;  % Vector from i1 to i3
        b32 = P3 - P2;  % Vector from i2 to i3
        b21 = P2 - P1;  % Vector from i1 to i2

        % Cross product to calculate orientation
        cc1 = cross(b13, b32);
        cc2 = cross(b21, b13);
        cc3 = cross(b32, b21);
        sun = (cc1 + cc2 + cc3)';  % Sum of cross products for orientation

        si = sign(sun * sv); % Calculate the sign for the triangle
        closure_value = (median_diff21 + median_diff32 - median_diff13) * si; % Calculate Closure Value

        % fprintf('   Triangle (i=%d, j=%d, k=%d):\n', i, j, k);
        % fprintf('     Baseline values: %.2f, %.2f, %.2f\n', median_diff21, median_diff32, median_diff13);
        % fprintf('     Initial Closure Value (with sign) = %.2f\n', closure_value);

        if abs(closure_value) < threshold % If the Closure value is within the threshold, no adjustment is needed
            % fprintf('     No adjustment needed. Closure value is already within the threshold.\n');
            continue;
        end

        adjustment_made = false;  % Flag to check if any adjustment happens

        for adjustment = [ambspace, -ambspace, 2*ambspace, -2*ambspace] % Try adjustments
            for baseline_idx = 1:3 % Attempt to correct each baseline by adjustments and check the new Closure value
                switch baseline_idx % Determine which baseline to adjust and the associated baseline ID
                    case 1
                        unfixedValue = median_diff21;
                        updateIdx = [i, j];
                    case 2
                        unfixedValue = median_diff32;
                        updateIdx = [j, k];
                    case 3
                        unfixedValue = median_diff13;
                        updateIdx = [i, k];
                end

                % Check current total adjustment in adjustmentsTable
                currentAdjustment = adjustmentsTable{updateIdx(1), updateIdx(2)};
                proposedAdjustment = currentAdjustment + adjustment;

                % Ensure the total adjustment remains within +-2*ambspace
                if proposedAdjustment > 2*ambspace || proposedAdjustment < -2*ambspace
                    continue;
                end

                tempBaseline = unfixedValue + adjustment; % Calculate the new baseline value after adjustment

                % Calculate the new Closure value with the adjusted baseline
                switch baseline_idx
                    case 1
                        new_closure_value = (tempBaseline + median_diff32 - median_diff13) * si;
                    case 2
                        new_closure_value = (median_diff21 + tempBaseline - median_diff13) * si;
                    case 3
                        new_closure_value = (median_diff21 + median_diff32 - tempBaseline) * si;
                end

                % Check if the new Closure value is within the threshold
                if abs(new_closure_value) < threshold
                    allBaselines(updateIdx(1), updateIdx(2)) = tempBaseline; % Update allBaselines with the new adjusted value
                    adjustmentsTable{updateIdx(1), updateIdx(2)} = proposedAdjustment; % Add the adjustment to the adjustmentsTable

                    % fprintf('     Adjustment saved in adjustmentsTable(%d, %d): %.2f\n', ...
                    %         updateIdx(1), updateIdx(2), adjustmentsTable{updateIdx(1), updateIdx(2)});

                    iteration_changes = true; % Mark a change
                    adjustment_made = true;
                    % fprintf('     Adjusted baseline (%d-%d): %.2f -> %.2f\n', updateIdx(1), updateIdx(2), unfixedValue, tempBaseline);
                    % fprintf('     New Closure Value = %.2f (Threshold satisfied)\n', new_closure_value);
                    break;  % Break once an adjustment is made
                end
            end
            if adjustment_made % If an adjustment was made, exit the outer loop as well
                break;
            end
        end

        if ~adjustment_made
            % fprintf('     No adjustment made for this baseline. Closure value not within threshold.\n');
        end
    end

    % if ~iteration_changes % If no changes were made, stop the iterations
    %     % fprintf('   No changes in this iteration. Stopping.\n');
    %     % break;
    % end
    
    if total_iteration_count == 5 % Break if there is endless loop problem with calculation
        break;
    end
    
end



%% Check if all closures are below the limit of half ambiguity spacing

closureLimitCheck = 0;

for idx = 1:size(combinations, 1)
    i = combinations(idx, 1);  % Station 1 (i1)
    j = combinations(idx, 2);  % Station 2 (i2)
    k = combinations(idx, 3);  % Station 3 (i3)
    updateIdx21 = [i, j];
    updateIdx32 = [j, k];
    updateIdx13 = [i, k];
    currentAdjustment21 = adjustmentsTable{updateIdx21(1), updateIdx21(2)};
    currentAdjustment32 = adjustmentsTable{updateIdx32(1), updateIdx32(2)};
    currentAdjustment13 = adjustmentsTable{updateIdx13(1), updateIdx13(2)};

    % Compute the Closure value only if there are no NaN values
    median_diff21 = allBaselines(i, j) + currentAdjustment21;  % i1 - i2
    median_diff32 = allBaselines(j, k) + currentAdjustment32;  % i2 - i3
    median_diff13 = allBaselines(i, k) + currentAdjustment13;  % i1 - i3

    % Check for NaN values and skip calculation if any are present
    if isnan(median_diff21) || isnan(median_diff32) || isnan(median_diff13)
        continue;  % Skip this combination if any value is NaN
    end

    % Centroid calculation
    P1 = [antenna(i).x, antenna(i).y, antenna(i).z]';
    P2 = [antenna(j).x, antenna(j).y, antenna(j).z]';
    P3 = [antenna(k).x, antenna(k).y, antenna(k).z]';
    sv = (P1 + P2 + P3) / 3;  % Centroid

    % Basis vectors from the stations
    b13 = P3 - P1;  % Vector from i1 to i3
    b32 = P3 - P2;  % Vector from i2 to i3
    b21 = P2 - P1;  % Vector from i1 to i2

    % Cross product to calculate orientation
    cc1 = cross(b13, b32);
    cc2 = cross(b21, b13);
    cc3 = cross(b32, b21);
    sun = (cc1 + cc2 + cc3)';  % Sum of cross products for orientation

    si = sign(sun * sv); % Calculate the sign for the triangle
    closure_value = (median_diff21 + median_diff32 - median_diff13) * si; % Calculate Closure Value

    if abs(closure_value) > 25 % If the Closure value is within the threshold, no adjustment is needed
        closureLimitCheck = closureLimitCheck+1;
    end
end

threshold = round(threshold);

if closureLimitCheck > 0
    fprintf('\t \t %d out of %d closure combinations are exceeding the closure limit of %d ns. \n', closureLimitCheck, length(combinations), threshold);
elseif closureLimitCheck == 0
    fprintf('\t \t All closures are below the closure limit of %d ns. \n',threshold);
end



%% Save Ambiguities into nscan and individual correction after automatic triangle delay closures
% nscan structure only for testing

nscan = scan;  % Copy of the original structure
dd = 1;

for row = 1:mnn  % If obs is saved as 4-3 in scan file and not 3-4
    for col = row+1:mnn % Only upper triangular part (row < col)
        % Mirror value from upper into lower triangle
        adjustmentsTable{col, row} = -adjustmentsTable{row, col};
        allBaselines(col, row) = -allBaselines(row, col);
    end
end

% Iterate through all scans and observations
for i = 1:length(scan)
    for j = 1:length(scan(i).obs)

        i1 = scan(i).obs(j).i1;
        i2 = scan(i).obs(j).i2;
        
        % Check if indices are within bounds
        if i1 <= size(adjustmentsTable, 1) && i2 <= size(adjustmentsTable, 2)
            adjustment_value = adjustmentsTable{i1, i2};
            if isnan(adjustment_value)
%                 fprintf('Skipping observation (%d, %d): Adjustment value is NaN.\n', i1, i2);
                continue;  % Skip if no ambiguity exists
            end

            % Retrieve the median value from allBaselines
            median_value = allBaselines(i1, i2);
            if isnan(median_value)
%                 fprintf('Skipping observation (%d, %d): Median value is NaN.\n', i1, i2);
                continue;  % Skip if no median value exists
            end
            
            % Compute the cases, important for individual correction after
            % automatic solution
            cases = [0, ambspace, -ambspace, 2*ambspace, -2*ambspace, 3*ambspace, -3*ambspace, 4*ambspace, -4*ambspace, 5*ambspace, -5*ambspace, 6*ambspace, -6*ambspace, 7*ambspace, -7*ambspace, 8*ambspace, -8*ambspace];
            cases2 = diffs(dd)+ adjustment_value + [0, ambspace, -ambspace, 2*ambspace, -2*ambspace, 3*ambspace, -3*ambspace, 4*ambspace, -4*ambspace, 5*ambspace, -5*ambspace, 6*ambspace, -6*ambspace, 7*ambspace, -7*ambspace, 8*ambspace, -8*ambspace];
            % cases = [0, ambspace, -ambspace, 2*ambspace, -2*ambspace];
            % cases2 = diffs(dd)+ adjustment_value + [0, ambspace, -ambspace, 2*ambspace, -2*ambspace];
         

            % Find the case closest to the median
            [~, closest_idx] = min(abs(cases2 - median_value));
            chosen_value = cases(closest_idx)+adjustment_value;
            
            % Calculate the final adjusted value
            final_adjustment = chosen_value;
%             final_adjustment = adjustment_value;
            
   
%             if abs(cases(closest_idx)) > 4*ambspace
%                 final_adjustment = adjustment_value;
%             end

            % Save the adjusted value into nscan
            if abs(final_adjustment) <= 7*abs(ambspace) && abs(diffs(dd)) < 1e+03
                nscan(i).obs(j).amb = final_adjustment * 10^-9;
                scan(i).obs(j).amb = final_adjustment * 10^-9; % Important so that nscan structure is not needed anymore
            end              
            dd = dd+1;
        end
    end
end

% Save the updated structure to file
% save([session, '_n_scan.mat'], 'nscan');
% fprintf('Adjusted values saved to %s_n_scan.mat\n', session);



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
        fprintf('\t \t Ambiguity list not available: %s\n', output_file_path);
    end
else
    fprintf('\t \t Ambiguities will not be changed \n');
end
% fprintf('\n')

file_id = fopen(output_file_path,'w+'); % Open/create new file for reading/writing, discard existing content

ambNumberCounter = 0;

% Iterate through all rows and columns of the adjustmentsTable
for i = 1:size(adjustmentsTable, 1)
    for j = 1:size(adjustmentsTable, 2)
        
        if isnan(adjustmentsTable{i, j}) % Skip NaN values 
            continue;
        end
        
        % Extract station names
        station1_name = strtrim(antenna(i).name); 
        station2_name = strtrim(antenna(j).name);
        
        % Iterate through all scans and observations to find matching i1 and i2
        for scan_idx = 1:length(scan)
            for obs_idx = 1:length(scan(scan_idx).obs)
                if scan(scan_idx).obs(obs_idx).i1 == i && ...
                   scan(scan_idx).obs(obs_idx).i2 == j
                    mjd = scan(scan_idx).mjd;
                    iso_idx = scan(scan_idx).iso; % ISO index for source name
                    source_name = strtrim(sources(iso_idx).name); % No extra spaces
                    % amb_value = nscan(scan_idx).obs(obs_idx).amb * 1e9; 
                    amb_value = scan(scan_idx).obs(obs_idx).amb * 1e9; % Use in order to be independet of nscan structure

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
                             
                    ambNumberCounter = ambNumberCounter + 1;         
                    
                    % Write the line to the file, using spaces for separation
                    fprintf(file_id, '%s\n', strjoin(line_data, ' ')); 
                end
            end
        end
    end
end

fclose(file_id);
fprintf('\t \t A total of %d ambiguities were found. \n', ambNumberCounter);
fprintf('\t \t Ambiguities written to %s\n', output_file_path);

% Count the number of NaN values and 0
% numNaNs = sum(isnan(diffs));
% numZeros = sum(diffs == 0);
% disp(['Number of NaN values: ', num2str(numNaNs)]);
% disp(['Number of zero values: ', num2str(numZeros)]);



%% Display Results
% Needed for testing only

% % Load data
% fileID = fopen(output_file_path, 'r'); 
% 
% if fileID == -1
%     error('File could not be opened.');
% end
% 
% % Read file data
% data = textscan(fileID, '%s %s %f %s %d', 'Delimiter', ' ', 'MultipleDelimsAsOne', 1);
% fclose(fileID);
% 
% % Extract ambiguity values
% ambiguity_values = data{5};
% 
% if isempty(ambiguity_values)
%     disp('No ambiguity values found.');
%     return;
% end
% 
% % Define the specific categories
% categories = [0, -50, -100, 50, 100];
% counts = zeros(size(categories));
% 
% % Count occurrences for defined categories
% for i = 1:length(categories)
%     counts(i) = sum(ambiguity_values == categories(i));
% end
% 
% % Count values outside defined categories
% counts_less_100 = sum(ambiguity_values < -100);
% counts_greater_100 = sum(ambiguity_values > 100);
% 
% % Create final table
% final_values = [categories, -101, 101]; % -101 and 101 are placeholders for < -100 and > 100
% final_counts = [counts, counts_less_100, counts_greater_100];

% disp('Adjustments Table (displaying adjustments made for baselines):');
% aad = round(adjustmentsTable);
% disp(aad);

% Display result
% ambiguity_table = table(final_values(:), final_counts(:), 'VariableNames', {'AmbNs', 'Count'});
% disp(ambiguity_table);

% valueCount0 = length(diffs) - length(ambiguity_values)
% 
% allsum = length(diffs)



%% Plot of Baselines
% Needed for testing only

% % tic;
% 
% % Unique (i1, i2) pairs
% uniquePairs = unique([[mb_sb_X.i1]', [mb_sb_X.i2]'], 'rows'); 
% 
% % Precompute indices and store
% i1Indices = [mb_sb_X.i1]';
% i2Indices = [mb_sb_X.i2]';
% diffs = [mb_sb_X.diff]';
% mjds = [mb_sb_X.mjd]';
% 
% % Analysis of each (i1, i2) pair
% for k = 1:size(uniquePairs, 1)
%     i1 = uniquePairs(k, 1);
%     i2 = uniquePairs(k, 2);
% 
%     % Extract entries for this (i1, i2) pair
%     indices = find(i1Indices == i1 & i2Indices == i2);
% 
%     % Skip if a baseline has less than 10 values
% %     if length(indices) < 10
% %         continue;
% %     end 
% 
%     % Extract values for diffs and mjds
%     diffsForPair = diffs(indices);
%     mjdsForPair = mjds(indices);
% 
%     % Create a separate plot for this (i1, i2) pair
%     if i1 == 4
%     % if i1 == 13 && i2 == 14
% 
%     figure;
%     plot(mjdsForPair, diffsForPair, '-o', 'LineWidth', 1.5, 'DisplayName', 'Data');
%     hold on;
%     % Initialize variables for legend entries
%     legendEntries = {'Data'};  % Start with the data plot in the legend
%     % Add grid and legend
%     grid on;
%     legend(legendEntries, 'Location', 'best');
%     title(sprintf('Baseline Differences: i1 = %d, i2 = %d', i1, i2));
%     xlabel('MJD (Modified Julian Date)');
%     ylabel('MB-SB difference (ns)');
%     hold off;
% 
%     % hold on;
% 
%     % figure;
%     % plot(mjdsForPair, sigObsDIFFForPair, '-o', 'LineWidth', 1.5, 'DisplayName', 'Data');
%     % hold on;
%     % % Initialize variables for legend entries
%     % legendEntries = {'Data'};  % Start with the data plot in the legend
%     % % Add grid and legend
%     % grid on;
%     % legend(legendEntries, 'Location', 'best');
%     % title(sprintf('Baseline Differences Sigma: i1 = %d, i2 = %d', i1, i2));
%     % xlabel('MJD (Modified Julian Date)');
%     % ylabel('MB-SB difference (ns)');
%     % hold off;
% 
% 
%     end
% 
% end
% 
% % toc;



%% Function to calculate the adjusted median

if medianSpecial == 1
    fprintf('\t \t WARNING - special case of median-calculation used, please check the observations\n'); 
end

% toc
endTime = round(toc);

fprintf('\t Ambig corr: finished after %d seconds \n', endTime')

function adjustedMedian = calculateAdjustedMedian(values, Amb)
    sortedValues = sort(values); % Sorted list of all values
    n = numel(sortedValues);

    % Median calculation based on odd or even number of observations
    if mod(n, 2) == 1 % Case for odd number
        medianValue = sortedValues((n + 1) / 2); % No further adjustment needed
    else % Case for even number
        medianValue = (sortedValues(n / 2) + sortedValues(n / 2 + 1)) / 2;
        first_value = sortedValues(n / 2);
        second_value = sortedValues(n / 2 + 1);
        valuediff = first_value - second_value;
        % Check if the difference is too big, special case (jump) where exactly
        % 50% of observations are shifted by ambiguities and 50 % are not
        if abs(valuediff) > 0.6*abs(Amb) % Check special case
            medianValue = sortedValues(n / 2);
            medianSpecial = 1;
        end   
    end
    
    adjustedMedian = medianValue;
end


end