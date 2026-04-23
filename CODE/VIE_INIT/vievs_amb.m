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
%       23.04.2026 by Peter Urban: major update of the approach for mixed
%       mode sessions, several bug fixes for the iterations and matrices
%
% ************************************************************************

function [amb1] = vievs_amb(out_struct,wrapper_data, parameter, antenna, sources, scan, fband_amb)


%% Load and prepare variables

tic;
format longg
fprintf('\t Ambig corr: directly calculated in VieVS \n')
amb1 = 1;
band = fband_amb; % choose the band for the calculation

if band == 'X'
    ambspaceV = (out_struct.Observables.AmbigSize_bX.AmbigSize.val).* 10^9; % Ambiguity spacing [ns]
elseif band == 'S'
    ambspaceV = (out_struct.Observables.AmbigSize_bS.AmbigSize.val).* 10^9; % Ambiguity spacing [ns]    
end

len_ambspaceV = length(ambspaceV);

ambspace = mode(ambspaceV);
if round(ambspace) == 250
    ambspace = ambspace/5;
end

threshold = ambspace/2;  % Threshold for closure value
maxIterationsAll = 30;
medianSpecial = 0; % check if a special case for median calculation exist
totalObservations=out_struct.head.NumObs.val; % Number of observations

% Multiband and Singleband Group Delay for X-Band
if band == 'X'
    SBD_X = double(out_struct.Observables.SBDelay_bX.SBDelay.val);
    SBD_X_sigma = double(out_struct.Observables.SBDelay_bX.SBDelaySig.val);
    MBD_X = double(out_struct.Observables.GroupDelay_bX.GroupDelay.val);
elseif band == 'S'
    SBD_X = double(out_struct.Observables.SBDelay_bS.SBDelay.val);
    MBD_X = double(out_struct.Observables.GroupDelay_bS.GroupDelay.val);   
    SBD_X_sigma = double(out_struct.Observables.SBDelay_bS.SBDelaySig.val);
end

diffX = (MBD_X - SBD_X).* 10^9; % mb-sb differences [ns]

idx_X = 1;
mb_sb_X(totalObservations) = struct('scan', [], 'obs', [], 'diff', [], 'ambSpace', [], 'i1', [], 'i2', [], 'mjd', [], 'isoV', [], 'qCodeX', [], 'qCodeS', [], 'SBsig', []);
mnn = length(antenna); % Total number of stations

% Read all the data and save it for calculations
for i = 1:length(scan) 
    for j = 1:length(scan(i).obs) 
        mb_sb_X(idx_X).scan = i;
        mb_sb_X(idx_X).obs = j;
        mb_sb_X(idx_X).diff = diffX(idx_X); % [ns]
        mb_sb_X(idx_X).SBsig = SBD_X_sigma(idx_X); % [s]
        if len_ambspaceV == 1
            mb_sb_X(idx_X).ambSpace = ambspaceV; % [ns]
        else
            mb_sb_X(idx_X).ambSpace = ambspaceV(idx_X); % [ns]
        end
        mb_sb_X(idx_X).i1 = scan(i).obs(j).i1;
        mb_sb_X(idx_X).i2 = scan(i).obs(j).i2;
        mb_sb_X(idx_X).qCodeX = str2double(scan(i).obs(j).q_code_X);
        mb_sb_X(idx_X).qCodeS = str2double(scan(i).obs(j).q_code_S);
        mb_sb_X(idx_X).i2 = scan(i).obs(j).i2;
        mb_sb_X(idx_X).mjd = scan(i).mjd;
        mb_sb_X(idx_X).isoV = scan(i).iso;
        idx_X = idx_X + 1;
    end
end



%% Ambiguity Calculation with Iterative Baseline Correction 

allBaselines = NaN(mnn, mnn);  % Table for all baseline medians
allAmbSpace = NaN(mnn, mnn);  % Table for all ambiguity spacing
allqCode = NaN(mnn, mnn); % Table for all quality codes
combinations = nchoosek(1:mnn, 3); % Generate all combinations of 3 stations 
i1Indices = [mb_sb_X.i1]';
i2Indices = [mb_sb_X.i2]';
diffs = [mb_sb_X.diff]';
ambSpaces = [mb_sb_X.ambSpace]';
qCX = [mb_sb_X.qCodeX]';
qCS = [mb_sb_X.qCodeS]';
SBsig = [mb_sb_X.SBsig]';

for i = 1:mnn
    for j = 1:mnn
        if i == j, continue; end  % Skip diagonal entries

        % Find all data for baseline i-j
        indices = find(i1Indices == i & i2Indices == j);
        if isempty(indices)
            continue;  % No data for this baseline
        end

        MBD_X_ns = MBD_X.* 10^9; % [ns]
        SBD_X_ns = SBD_X.* 10^9; % [ns]
        MBD_X_pair_ns = MBD_X_ns(indices); % [ns]
        SBD_X_pair_ns = SBD_X_ns(indices); % [ns]
        diffsForPair = diffs(indices);
        ambSpacePair = ambSpaces(indices);
        qCXpair = qCX(indices);
        SbSigpair = SBsig(indices);
        qCSpair = qCS(indices);

        % Calculate median values for the baselines
        [adjustedMedian, qCodeMedian] = calculateAdjustedMedian(i, j, antenna, diffsForPair, qCXpair, ambspace, SbSigpair, MBD_X_pair_ns, SBD_X_pair_ns);

        adjustedMedianAmb = median(ambSpacePair);
        allBaselines(i, j) = adjustedMedian;
        allAmbSpace(i, j) = adjustedMedianAmb;
        allqCode(i, j) = qCodeMedian;

        if ~isnan(adjustedMedian) && isnan(adjustedMedianAmb)
            allAmbSpace(i, j) = ambspace;
        end
    end
end

% Table for the storage of ambiguity values
adjustmentsTable = array2table(allBaselines, 'VariableNames', arrayfun(@(x) sprintf('S%d', x), 1:mnn, 'UniformOutput', false));

% Replace NaN values in adjustmentsTable with NaN and all non-NaN values with 0
for row = 1:mnn
    for col = 1:mnn
        if isnan(allBaselines(row, col)) || allBaselines(row, col) == 0
            if row < col && ~isnan(allBaselines(col, row))
                % If in the upper triangle (row < col) NaN, but there is a value in lower triangle (col, row) 
                allBaselines(row, col) = allBaselines(col, row);
                allqCode(row, col) = allqCode(col, row);
                allAmbSpace(row, col) = allAmbSpace(col, row);
                adjustmentsTable{col, row} = 0;
            end
        end
        adjustmentsTable{row, col} = 0;
        adjustmentsTable{row, col} = 0;
    end
end

for row = 2:mnn % nan values for lower left triangle
    for col = 1:row-1
        adjustmentsTable{row, col} = NaN;
        allBaselines(row, col) = NaN;
        allqCode(row, col) = NaN;
    end
end

ALLBASELINES_start = round(allBaselines);
Check_MATRIX = ones(size(ALLBASELINES_start));
iteration_count = 0;  % Counter for number of iterations
allBaselines_best = allBaselines;
adjustmentsTable_best = adjustmentsTable;
closureCounter_best = 10000;
allBaselines_new = allBaselines;
adjustmentsTable_new = adjustmentsTable;

% Main calculation of ambiguities with closures
for ixc = 1:maxIterationsAll % calculate 30 solutions with random order of the triangles (pick best one later)
    total_iteration_count = 0;
    allBaselines = allBaselines_new;
    adjustmentsTable = adjustmentsTable_new;
    
    while total_iteration_count < 10 % 10x calculation of whole network within one solution
        iteration_count = iteration_count + 1;
        total_iteration_count = total_iteration_count + 1;
    
        if ixc > 1
            combinations = nchoosek(1:mnn, 3); % Generate all combinations of 3 stations 
            comb_rand = combinations(randperm(size(combinations,1)), :);
            combinations = comb_rand;
        end

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
                continue;  % Skip this combination if any value is NaN
            end

            Check_MATRIX_21 = Check_MATRIX(i, j);
            Check_MATRIX_32 = Check_MATRIX(j, k);
            Check_MATRIX_13 = Check_MATRIX(i, k);

            if Check_MATRIX_21 == 0
                continue;
            elseif Check_MATRIX_32 == 0
                continue;
            elseif Check_MATRIX_13 == 0
                continue;
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

            threshold = ( (allAmbSpace(i, j)+allAmbSpace(j, k)+allAmbSpace(i, k) ) / 3 ) / 2;

            if abs(closure_value) < threshold % If the Closure value is within the threshold, no adjustment is needed
                continue;
            end

            iteration_changes = false;

            ambspace1 = allAmbSpace(i, j); 
            % Try ambiguity correction for first baseline
            for adjustment = [ambspace1, -ambspace1, 2*ambspace1, -2*ambspace1, 3*ambspace1, -3*ambspace1, 4*ambspace1, -4*ambspace1, 5*ambspace1, -5*ambspace1] % Try adjustments
                unfixedValue = median_diff21;
                updateIdx = [i, j];
                currentAdjustment = adjustmentsTable{updateIdx(1), updateIdx(2)};
                proposedAdjustment = currentAdjustment + adjustment;
                tempBaseline = unfixedValue + adjustment; 
                new_closure_value = (tempBaseline + median_diff32 - median_diff13) * si;
                if Check_MATRIX(i,j) == 1 && Check_MATRIX(j,k) == 1 && Check_MATRIX(i,k) == 1
                    if abs(new_closure_value) < threshold  
                        allBaselines_test = allBaselines;

                        combinations1 = nchoosek(1:mnn, 3); % Generate all combinations of 3 stations 
                        comb_rand1 = combinations1(randperm(size(combinations1,1)), :);
                        combinations1 = comb_rand1;

                        [closureLimitCheck_norm1] = calcClosures(combinations1, allBaselines, antenna, threshold);
                        allBaselines_test(updateIdx(1), updateIdx(2)) = tempBaseline; % Update allBaselines with the new adjusted value
                        [closureLimitCheck_test1] = calcClosures(combinations1, allBaselines_test, antenna, threshold);
                
                        if closureLimitCheck_test1 < closureLimitCheck_norm1 && maxIterationsAll <= 15
                            allBaselines(updateIdx(1), updateIdx(2)) = tempBaseline; % Update allBaselines with the new adjusted value
                            adjustmentsTable{updateIdx(1), updateIdx(2)} = proposedAdjustment; % Add the adjustment to the adjustmentsTable
                            iteration_changes = true; % Mark a change
                            break;  % Break once an adjustment is made
                        elseif closureLimitCheck_test1 <= closureLimitCheck_norm1 && maxIterationsAll > 15
                            allBaselines(updateIdx(1), updateIdx(2)) = tempBaseline; % Update allBaselines with the new adjusted value
                            adjustmentsTable{updateIdx(1), updateIdx(2)} = proposedAdjustment; % Add the adjustment to the adjustmentsTable
                            iteration_changes = true; % Mark a change
                            break;  % Break once an adjustment is made
                        end
                    end
                end
            end

            if iteration_changes == false % Try ambiguity correction for second baseline
                ambspace2 = allAmbSpace(j, k);
                for adjustment = [ambspace2, -ambspace2, 2*ambspace2, -2*ambspace2, 3*ambspace2, -3*ambspace2, 4*ambspace2, -4*ambspace2, 5*ambspace2, -5*ambspace2] % Try adjustments
                    unfixedValue = median_diff32;
                    updateIdx = [j, k];
                    currentAdjustment = adjustmentsTable{updateIdx(1), updateIdx(2)};
                    proposedAdjustment = currentAdjustment + adjustment;
                    tempBaseline = unfixedValue + adjustment; 
                    new_closure_value = (median_diff21 + tempBaseline - median_diff13) * si;

                    if Check_MATRIX(i,j) == 1 && Check_MATRIX(j,k) == 1 && Check_MATRIX(i,k) == 1
                        if abs(new_closure_value) < threshold
                            allBaselines_test = allBaselines;

                            combinations2 = nchoosek(1:mnn, 3); % Generate all combinations of 3 stations 
                            comb_rand2 = combinations2(randperm(size(combinations2,1)), :);
                            combinations2 = comb_rand2;

                            [closureLimitCheck_norm2] = calcClosures(combinations2, allBaselines, antenna, threshold);
                            allBaselines_test(updateIdx(1), updateIdx(2)) = tempBaseline; % Update allBaselines with the new adjusted value
                            [closureLimitCheck_test2] = calcClosures(combinations2, allBaselines_test, antenna, threshold);

                            if closureLimitCheck_test2 < closureLimitCheck_norm2 && maxIterationsAll <= 15
                                allBaselines(updateIdx(1), updateIdx(2)) = tempBaseline; % Update allBaselines with the new adjusted value
                                adjustmentsTable{updateIdx(1), updateIdx(2)} = proposedAdjustment; % Add the adjustment to the adjustmentsTable
                                iteration_changes = true; % Mark a change
                                break;  % Break once an adjustment is made
                            elseif closureLimitCheck_test2 <= closureLimitCheck_norm2 && maxIterationsAll > 15
                                allBaselines(updateIdx(1), updateIdx(2)) = tempBaseline; % Update allBaselines with the new adjusted value
                                adjustmentsTable{updateIdx(1), updateIdx(2)} = proposedAdjustment; % Add the adjustment to the adjustmentsTable
                                iteration_changes = true; % Mark a change
                                break;  % Break once an adjustment is made
                            end
                        end 
                    end
                end
            end

            if iteration_changes == false % Try ambiguity correction for third baseline
                ambspace3 = allAmbSpace(i, k);
                for adjustment = [ambspace3, -ambspace3, 2*ambspace3, -2*ambspace3, 3*ambspace3, -3*ambspace3, 4*ambspace3, -4*ambspace3, 5*ambspace3, -5*ambspace3] % Try adjustments
                    unfixedValue = median_diff13;
                    updateIdx = [i, k];
                    currentAdjustment = adjustmentsTable{updateIdx(1), updateIdx(2)};
                    proposedAdjustment = currentAdjustment + adjustment;
                    tempBaseline = unfixedValue + adjustment; 
                    new_closure_value = (median_diff21 + median_diff32 - tempBaseline) * si;

                    if Check_MATRIX(i,j) == 1 && Check_MATRIX(j,k) == 1 && Check_MATRIX(i,k) == 1
                        if abs(new_closure_value) < threshold
                            allBaselines_test = allBaselines;

                            combinations3 = nchoosek(1:mnn, 3); % Generate all combinations of 3 stations 
                            comb_rand3 = combinations3(randperm(size(combinations3,1)), :);
                            combinations3 = comb_rand3;

                            [closureLimitCheck_norm3] = calcClosures(combinations3, allBaselines, antenna, threshold);
                            allBaselines_test(updateIdx(1), updateIdx(2)) = tempBaseline; % Update allBaselines with the new adjusted value
                            [closureLimitCheck_test3] = calcClosures(combinations3, allBaselines_test, antenna, threshold);
                    
                            if closureLimitCheck_test3 < closureLimitCheck_norm3 && maxIterationsAll <= 15
                                allBaselines(updateIdx(1), updateIdx(2)) = tempBaseline; % Update allBaselines with the new adjusted value
                                adjustmentsTable{updateIdx(1), updateIdx(2)} = proposedAdjustment; % Add the adjustment to the adjustmentsTable
                                iteration_changes = true; % Mark a change
                                break;  % Break once an adjustment is made
                            elseif closureLimitCheck_test3 <= closureLimitCheck_norm3 && maxIterationsAll > 15
                                allBaselines(updateIdx(1), updateIdx(2)) = tempBaseline; % Update allBaselines with the new adjusted value
                                adjustmentsTable{updateIdx(1), updateIdx(2)} = proposedAdjustment; % Add the adjustment to the adjustmentsTable
                                iteration_changes = true; % Mark a change
                                break;  % Break once an adjustment is made
                            end
                        end
                    end
                end
            end

        end
    
        closureLimitCheck4 = 0;
        [closureLimitCheck4] = calcClosures(combinations, allBaselines, antenna, threshold);

        if total_iteration_count == 10 % Break if there is endless loop problem with calculation
            break;
        end

    end

    if closureCounter_best > closureLimitCheck4
        allBaselines_best = allBaselines;
        adjustmentsTable_best = adjustmentsTable;
        closureCounter_best = closureLimitCheck4;
    end

    closureCounter(ixc) = closureLimitCheck4;
    if closureLimitCheck4 == 0
        break;
    end

end % end the first for loop

allBaselines = allBaselines_best;
adjustmentsTable = adjustmentsTable_best;

% Check if all closures are below the limit of half ambiguity spacing
threshold = round(threshold);
if closureCounter_best > 0
    fprintf('\t \t %d out of %d closure combinations are exceeding the closure limit of %d ns. \n', closureCounter_best, length(combinations), threshold);
elseif closureCounter_best == 0
    fprintf('\t \t All closures are below the closure limit of %d ns. \n',threshold);
end



%% Corrections within a baseline relative to the median diff value of the baseline

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
                continue;  % Skip if no ambiguity exists
            end

            % Retrieve the median value from allBaselines
            median_value = allBaselines(i1, i2);
            adAmbSpace = allAmbSpace(i1, i2);
            if isnan(median_value)
                continue;  % Skip if no median value exists
            end

            diffsA = [mb_sb_X.diff]';
            
            cases = [0, adAmbSpace, -adAmbSpace, 2*adAmbSpace, -2*adAmbSpace, 3*adAmbSpace, -3*adAmbSpace, 4*adAmbSpace, -4*adAmbSpace];
            cases2 = diffsA(dd)+ adjustment_value + [0, adAmbSpace, -adAmbSpace, 2*adAmbSpace, -2*adAmbSpace, 3*adAmbSpace, -3*adAmbSpace, 4*adAmbSpace, -4*adAmbSpace];

            % Find the case closest to the median
            [~, closest_idx] = min(abs(cases2 - median_value));
            chosen_value = cases(closest_idx)+adjustment_value;
            
            % Calculate the final adjusted value
            final_adjustment = chosen_value;

            % Save the adjusted value 
            if abs(final_adjustment) <= 7*abs(ambspace) && abs(diffsA(dd)) < 1e+03
                scan(i).obs(j).amb = final_adjustment * 10^-9; 
            end              
            dd = dd+1;
        end
    end
end



%% Save Ambiguities to TXT File

% Ambiguity file
parameter.amb.flag_change_amb = true; %true

checkPath = ['../DATA/AMB/', parameter.filepath(end-4:end-1), '/'];
if ~exist(checkPath,'dir')
    mkdir(checkPath);
end

output_file_path = ['../DATA/AMB/', parameter.filepath(end-4:end-1), '/', [parameter.session_name '_' parameter.vie_init.vgosDb_observation_parameter(end) ], '.AMB'];

exist(output_file_path, 'file');

if parameter.amb.flag_change_amb
    if exist(output_file_path, 'file')
    else
        fprintf('\t \t Ambiguity list not available: %s\n', output_file_path);
    end
else
    fprintf('\t \t Ambiguities will not be changed \n');
end
% fprintf('\n')

file_id = fopen(output_file_path,'w+'); % Open/create new file for reading/writing, discard existing content

ambNumberCounter = 0;

% Iterate through all scans and observations
for scan_idx = 1:length(scan)
    for obs_idx = 1:length(scan(scan_idx).obs)
        
        i = scan(scan_idx).obs(obs_idx).i1;
        j = scan(scan_idx).obs(obs_idx).i2;
        
        % Extract station names
        station1_name = strtrim(antenna(i).name); 
        station2_name = strtrim(antenna(j).name);
        
        mjd = scan(scan_idx).mjd;
        iso_idx = scan(scan_idx).iso;
        source_name = strtrim(sources(iso_idx).name);
        
        amb_value = scan(scan_idx).obs(obs_idx).amb * 1e9;

        formatted_mjd = sprintf('%.12f', mjd);
        
        if amb_value > 0
            formatted_corr = sprintf('%d', round(amb_value));
            ambNumberCounter = ambNumberCounter + 1;
        elseif amb_value == 0
            formatted_corr = sprintf('%d', round(amb_value));
        else
            formatted_corr = sprintf('-%d', abs(round(amb_value)));
            ambNumberCounter = ambNumberCounter + 1;
        end

        % Format output line
        line_data = {sprintf('%-8s', station1_name), ...
                     sprintf('%-8s', station2_name), ...
                     sprintf('%-18s', formatted_mjd), ...
                     sprintf('%-8s', source_name), ...
                     sprintf('%-5s', formatted_corr)};

        fprintf(file_id, '%s\n', strjoin(line_data, ' ')); 
    end
end

fclose(file_id);
fprintf('\t \t A total of %d ambiguities were found. \n', ambNumberCounter);
fprintf('\t \t Ambiguities written to %s\n', output_file_path);



%% Function to calculate the adjusted median and closures

if medianSpecial == 1
    fprintf('\t \t WARNING - special case of median-calculation used, please check the observations\n'); 
end

endTime = round(toc);
fprintf('\t Ambig corr: finished after %.1f seconds \n', endTime')



function [adjustedMedian, qCodeMedian] = calculateAdjustedMedian(i1, i2, antenna, values, qCXpair, Amb, SbSigpair, MBD_X_pair_ns, SBD_X_pair_ns)
    
    id9 = find(qCXpair == 9 & SbSigpair < 1e-8);
    id8 = find(qCXpair == 8 & SbSigpair < 1e-8);
    id7 = find(qCXpair == 7 & SbSigpair < 1e-8);
    id6 = find(qCXpair == 6 & SbSigpair < 1e-8);
    id5 = find(qCXpair == 5 & SbSigpair < 1e-8);
    id4 = find(qCXpair == 4 & SbSigpair < 1e-8);
    id3 = find(qCXpair == 3 & SbSigpair < 1e-8);
    id2 = find(qCXpair == 2 & SbSigpair < 1e-8);
    id1 = find(qCXpair == 1 & SbSigpair < 1e-8);
    id0 = find(qCXpair == 0 & SbSigpair < 1e-8);
    if ~isempty(id9) 
        idd = id9;
    elseif ~isempty(id8) 
        idd = id8;
    elseif ~isempty(id7) 
        idd = id7;
    elseif ~isempty(id6) 
        idd = id6;
    elseif ~isempty(id5) 
        idd = id5;
    elseif ~isempty(id4) 
        idd = id4;
    elseif ~isempty(id3) 
        idd = id3;
    elseif ~isempty(id2) 
        idd = id2;
    elseif ~isempty(id1) 
        idd = id1;
    elseif ~isempty(id0) 
        idd = id0;
    else % if there is a general problem
        idd = 111;
    end 
    
    if idd ~= 111
        values = values(idd);
        qCXpair = qCXpair(idd);
    end

    sortedValues = sort(values); % Sorted list of all values
    n = numel(sortedValues);

    % Median calculation based on odd or even number of observations
    if mod(n, 2) == 1 % Case for odd number
        medianValue = sortedValues((n + 1) / 2); % No further adjustment needed
        idxxx = find(values == medianValue);
        qCodeMedian = qCXpair(idxxx);
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
        idxxx1 = find(values == first_value);
        idxxx2 = find(values == second_value);
        if idxxx1 < idxxx2
            qCodeMedian = qCXpair(idxxx1);
        elseif medianValue == 0
            qCodeMedian = 0;
        else
            qCodeMedian = qCXpair(idxxx2);
        end
    end

    qCodeMedian = qCodeMedian(1);
    adjustedMedian = medianValue;

end



function [closureLimitCheck] = calcClosures(combinations, allBaselines, antenna, threshold)
    closureLimitCheck = 0;
    for idxxxx = 1:size(combinations, 1)
        i = combinations(idxxxx, 1);  % Station 1 (i1)
        j = combinations(idxxxx, 2);  % Station 2 (i2)
        k = combinations(idxxxx, 3);  % Station 3 (i3)

        median_diff21 = allBaselines(i, j);  % i1 - i2
        median_diff32 = allBaselines(j, k);  % i2 - i3
        median_diff13 = allBaselines(i, k);  % i1 - i3

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

        if abs(closure_value) > threshold % If the Closure value is within the threshold, no adjustment is needed
            closureLimitCheck = closureLimitCheck+1;
        end
    end
end



end







