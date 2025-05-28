% ************************************************************************
%   Description:
%   This function calculate the ionospheric contribution and the
%   ionospheric flag.
%
%   Input:	
%      Variables from out_struct, wrapper_data
%
%   Output:
%      Calculated variables: ionospheric contribution, 
%      ionospheric formal error and quality-flag
%
%   Coded for VieVS: 
%   21 Nov 2024 by Peter Urban
%
%   Update:
%   07 May 2025 by Peter Urban: sign for hBw got switched, usb and lsb
%   lines for data reading got switched, in total results dont change to
%   previous version
%   19 May 2025 by Peter Urban: Number of bits per sample added, correction
%   of various errors
%   21 May 2025 by Peter Urban: Number of bits per sample removed - it 
%   should always be samplerate/4, only positive samplerate for
%   calculation, samplerate special cases 
%   28 May 2025 by Peter Urban: Small Updates, iono-flag if QCode S/X = 0,
%   channel weight stays 0 and will not be set to 1, channel frequency
%   check for repeat attribute if there are not values for every
%   observation, set NaN values to 0 for effective frequency and iono
%   corretion
%
% ************************************************************************

function [iono_x_corr_own, sigma_iono_x_corr_own, qflag_ion_own] = vievs_iono(out_struct,wrapper_data)

% 1) Define all the needed variables

try % Try to get length for obs in oder to have length for flags if there is a problem
    GroupDelayX = out_struct.Observables.GroupDelay_bX.GroupDelay.val;
    lenOBS = length(GroupDelayX);

catch

end

check = 0; % variable for checking if all needed data was read correct

try % load all the variables from the vgosdb
    QualityCode_Xv = out_struct.Observables.QualityCode_bX.QualityCode.val;
    QualityCode_X = double(QualityCode_Xv) - 48; % Quality Code Flag
    ChanAmpPhase_X = out_struct.Observables.ChannelInfo_bX.ChanAmpPhase.val; % amplitude and phase of the channels
    ChannelFreq_X = out_struct.Observables.ChannelInfo_bX.ChannelFreq.val; % channel frequency [MHz] 
    NumChannels_X = out_struct.Observables.ChannelInfo_bX.NumChannels.val; % number of channels
    NumSamples_X = out_struct.Observables.ChannelInfo_bX.NumSamples.val; % number of samples by sideband and channel 
    RefFreq_X = out_struct.Observables.RefFreq_bX.RefFreq.val; % reference frequency [MHz]
    SampleRate_X = out_struct.Observables.ChannelInfo_bX.SampleRate.val; % samplerate [Hz]
    SampleRate_X = SampleRate_X(1);
    if SampleRate_X < 0
        fprintf('WARNING vievs_iono: negative samplerate in vgosDB - positive value used for calculation \n')
    end
    if SampleRate_X == -32768000 % problem with samplerate in r1 ~ 2010
        SampleRate_X = 16000000;
    elseif SampleRate_X == 9216000 % problem with samplerate in r1 ~ 2010
        SampleRate_X = 16000000;
    elseif SampleRate_X == 64000000
        try  % Check for station NYALE13N - assume 16 MHZ samplerate
            nn=size(out_struct.head.StationList.val,2);
            for ii=1:nn
                curName=out_struct.head.StationList.val(:,ii)';
                curNameLong=[curName, '          '];
                stati(ii).name=curNameLong(1:8);
            end
            names = {stati.name};
            checkST = contains(names,'NYALE13N');
            scheckST = sum(checkST); 
            if scheckST == 1 % station NYALE13N included
                SampleRate_X = 16000000;
            end
        catch
        end
    end
    SampleRate_X = abs(SampleRate_X); % samplerate only has positive values
    % BitSample_X = out_struct.Observables.ChannelInfo_bX.BITSAMPL.val; % Number of bits per sample
    % BitSample_X = max(BitSample_X); % If there is saved one value for each observation
    HalfBwX = SampleRate_X/4/1.0e6; % half bandwidth [MHz]
    Repeat_X = out_struct.Observables.ChannelInfo_bX.ChannelFreq.attr; % check if channel frequency is same for all observations
    if length(ChannelFreq_X) < length(QualityCode_X) 
        if Repeat_X(3).name == 'REPEAT'
            % in this case the channel frequency stays the same for all observations
        else
            fprintf('WARNING vievs_iono: Channel Frequency not available for every obs - Repeat is missing in vgosDB - same Channel Frequency used for every obs \n')
        end
    end

    QualityCode_Sv = out_struct.Observables.QualityCode_bS.QualityCode.val;
    QualityCode_S = double(QualityCode_Sv) - 48; % Quality Code Flag
    ChanAmpPhase_S = out_struct.Observables.ChannelInfo_bS.ChanAmpPhase.val; % amplitude and phase of the channels
    ChannelFreq_S = out_struct.Observables.ChannelInfo_bS.ChannelFreq.val; % channel frequency [MHz] 
    NumChannels_S = out_struct.Observables.ChannelInfo_bS.NumChannels.val; % number of channels
    NumSamples_S = out_struct.Observables.ChannelInfo_bS.NumSamples.val; % number of samples by sideband and channel 
    RefFreq_S = out_struct.Observables.RefFreq_bS.RefFreq.val; % reference frequency [MHz]
    SampleRate_S = out_struct.Observables.ChannelInfo_bS.SampleRate.val; % samplerate [Hz]
    SampleRate_S = SampleRate_S(1);
    if SampleRate_S < 0
        fprintf('WARNING vievs_iono: negative samplerate in vgosDB - positive value used for calculation \n')
    end
    if SampleRate_S == -32768000 % problem with samplerate in r1 ~ 2010
        SampleRate_S = 16000000;
    elseif SampleRate_S == 9216000 % problem with samplerate in r1 ~ 2010
        SampleRate_S = 16000000;
    elseif SampleRate_S == 64000000
        try  % Check for station NYALE13N - assume 16 MHZ samplerate
            nn=size(out_struct.head.StationList.val,2);
            for ii=1:nn
                curName=out_struct.head.StationList.val(:,ii)';
                curNameLong=[curName, '          '];
                stati(ii).name=curNameLong(1:8);
            end
            names = {stati.name};
            checkST = contains(names,'NYALE13N');
            scheckST = sum(checkST); 
            if scheckST == 1 % station NYALE13N included
                SampleRate_S = 16000000;
            end
        catch
        end
    end
    SampleRate_S = abs(SampleRate_S); % samplerate only has positive values
    HalfBwS = SampleRate_S/4/1.0e6; % half bandwidth [MHz]
    Repeat_S = out_struct.Observables.ChannelInfo_bS.ChannelFreq.attr; % check if channel frequency is same for all observations
    if length(ChannelFreq_S) < length(QualityCode_S) 
        if Repeat_S(3).name == 'REPEAT'
            % in this case the channel frequency stays the same for all observations
        else
            fprintf('WARNING vievs_iono: Channel Frequency not available for every obs - Repeat is missing in vgosDB - same Channel Frequency used for every obs \n')
        end
    end
    
    % get the correct name for the GroupDelayFull in all different cases
    tau_xx = wrapper_data.Observation.ObsEdit.files; 
    tauxx = tau_xx{end};  % extract the last filename from the cell array
    parts = split(tauxx, {'_', '.'});  % split by underscore and dot
    identifier = parts{2};  % get the analysis center from the filename
    type = 'GroupDelayFull';  
    identifier = convertCharsToStrings(identifier);
    idx = "bX";
    ids = "bS";
     
    % get the variable with the correct analysis center
    if identifier == idx || identifier == ids
        variable_name_bS = [type, '_bS'];  % e.g., 'GroupDelayFull_iIVS_bS'
        variable_name_bX = [type, '_bX'];  % e.g., 'GroupDelayFull_iIVS_bX'
    else
        variable_name_bS = [type, '_', identifier, '_bS'];  % e.g., 'GroupDelayFull_iIVS_bS'
        variable_name_bX = [type, '_', identifier, '_bX'];  % e.g., 'GroupDelayFull_iIVS_bX'
        variable_name_bS = strjoin(variable_name_bS, '');
        variable_name_bX = strjoin(variable_name_bX, '');
    end
    
    % variable_name_bS = strjoin(variable_name_bS, '')
    % variable_name_bX = strjoin(variable_name_bX, '')
    
    % obs values where the ambiguities are already included
    tau_x = out_struct.ObsEdit.(variable_name_bX).GroupDelayFull.val; % [s] 
    tau_s = out_struct.ObsEdit.(variable_name_bS).GroupDelayFull.val; % [s] 

    % Delay Measurement Sigma (Group Delay)
    sigma_tau_x = out_struct.Observables.GroupDelay_bX.GroupDelaySig.val; % [s] 
    sigma_tau_s = out_struct.Observables.GroupDelay_bS.GroupDelaySig.val; % [s] 

catch
    check = 1;
    fprintf('WARNING vievs_iono: No Calculation due to lack of availability of all variables \n')
    iono_x_corr_own = 0;
    sigma_iono_x_corr_own = 0;
    qflag_ion_own = zeros(lenOBS, 1);
    qflag_ion_own(qflag_ion_own == 0) = -1;
end


if check == 0

    try    
        % 2) Effective Frequency for X- and S-Band

        % total number of observations
        num_obs = size(ChanAmpPhase_X, 3);
        
        % check if usb lsb got switched later after May 2018
        checkLines = NumSamples_S(2, 1:num_obs); % second line for usb
        scheckLines = sum(checkLines);
        checkUSBLSB = 0;
        if scheckLines == 0
            checkUSBLSB = 1;
        end

        % pre-allocate a result array
        vx = zeros(num_obs, 1);
        vs = zeros(num_obs, 1);

        % calculate formal errors for effective frequencies
        sigma_vx = zeros(num_obs, 1);
        sigma_vs = zeros(num_obs, 1);

        % quality flag for the ionospheric contribution
        qflag_ion_own = zeros(num_obs, 1);

        % effective frequency for X-Band
        for i = 1:num_obs 
            if length(NumChannels_X)>1
                n = double(NumChannels_X(i)); % number of channels 
            else
                n = double(NumChannels_X); % number of channels 
            end
            ri = ChanAmpPhase_X(1, 1:n, i); % channel amplitudes (first line of ChanAmpPhase)

            if checkUSBLSB == 0 
                usb = NumSamples_X(2, 1:n, i); % second line for usb
                lsb = NumSamples_X(1, 1:n, i); % first line for lsb
            end
            if checkUSBLSB == 1 % switch for sessions after May 2018
                usb = NumSamples_X(1, 1:n, i); % first line for usb
                lsb = NumSamples_X(2, 1:n, i); % second line for lsb
            end

            if length(ChannelFreq_X) > 20
                vi = ChannelFreq_X(1:n, i); % channel frequency [MHz]
            else
                vi = ChannelFreq_X; % channel frequency [MHz] 
            end
            if length(RefFreq_X) > 5
                v0 = sum(double(RefFreq_X))/length(RefFreq_X);
            else
                v0 = double(RefFreq_X); % reference frequency [MHz]
            end

            wi = ((usb + lsb) .* ri); % weight for sums

            % adjust frequencies for USB/LSB confusion
            for j = 1:n
                if usb(j) > 0.0 && lsb(j) > 0.0
                elseif usb(j) > 0.0 && lsb(j) == 0.0      
                    vi(j) = vi(j) + HalfBwX; % add half bandwith for usb
                elseif lsb(j) > 0.0 && usb(j) == 0.0       
                    vi(j) = vi(j) - HalfBwX; % subtract half bandwith for lsb
                end
            end

            % computation of the sums
            sum_ri = sum(wi);
            sum_ri_vivi = wi * (vi - v0).^2;
            sum_ri_vi_v0 = wi * (vi - v0);
            sum_ri_vi =  wi * (1./vi);
            sum_ri_vi_v0_vi = wi * ((vi - v0) ./ vi);
    
            numerator = sum_ri * sum_ri_vivi - (sum_ri_vi_v0)^2;
            denominator = sum_ri_vi_v0 * sum_ri_vi - sum_ri * sum_ri_vi_v0_vi;

            vx(i) = sqrt(numerator / denominator); % effective frequency X-Band [MHz]

            % calculate formal error for vx
            sigma_vx(i) = (std(wi) / sqrt(n))./10^6; % [MHz] 
            % formal errors should be < 0.69 [MHz]
    
            % iono-flag where vx could not be calculated in the correct way
            if isnan(vx(i)) || vx(i) == 0
                qflag_ion_own(i) = -1;
                vx(i) = 0;
            elseif QualityCode_X(i) == 0
                qflag_ion_own(i) = -1;
            end
        end

        % effective frequency for S-Band
        for i = 1:num_obs 
            if length(NumChannels_S)>1
                n = double(NumChannels_S(i)); % number of channels 
            else
                n = double(NumChannels_S); % number of channels 
            end
            ri = ChanAmpPhase_S(1, 1:n, i); % channel amplitudes (first line of ChanAmpPhase)

            if checkUSBLSB == 0 
                usb = NumSamples_S(2, 1:n, i); % second line for usb
                lsb = NumSamples_S(1, 1:n, i); % first line for lsb
            end
            if checkUSBLSB == 1 % switch for sessions after May 2018 
                usb = NumSamples_S(1, 1:n, i); % first line for usb
                lsb = NumSamples_S(2, 1:n, i); % second line for lsb
            end

            if length(ChannelFreq_S) > 20
                vi = ChannelFreq_S(1:n, i); % channel frequency [MHz]
            else
                vi = ChannelFreq_S; % channel frequency [MHz] 
            end
            if length(RefFreq_S) > 5
                v0 = sum(double(RefFreq_S))/length(RefFreq_S);
            else
                v0 = double(RefFreq_S); % reference frequency [MHz]
            end

            wi = ((usb + lsb) .* ri); % weight for sums

            % adjust frequencies for USB/LSB confusion
            for j = 1:n
                if usb(j) > 0.0 && lsb(j) > 0.0
                elseif usb(j) > 0.0 && lsb(j) == 0.0
                    vi(j) = vi(j) + HalfBwS; % add half bandwith for usb
                elseif lsb(j) > 0.0 && usb(j) == 0.0    
                    vi(j) = vi(j) - HalfBwS; % subtract half bandwith for lsb
                end
            end

            % computation of the sums
            sum_ri = sum(wi);
            sum_ri_vivi = wi * (vi - v0).^2;
            sum_ri_vi_v0 = wi * (vi - v0);
            sum_ri_vi =  wi * (1./vi);
            sum_ri_vi_v0_vi = wi * ((vi - v0) ./ vi);
    
            numerator = sum_ri * sum_ri_vivi - (sum_ri_vi_v0)^2;
            denominator = sum_ri_vi_v0 * sum_ri_vi - sum_ri * sum_ri_vi_v0_vi;

            vs(i) = sqrt(numerator / denominator); % effective frequency S-Band [MHz]

            % calculate formal error for vs
            sigma_vs(i) = (std(wi) / sqrt(n))./10^6; % [MHz]  
            % formal errors should be < 0.26 [MHz]
        
            % iono-flag where vs could not be calculated in the correct way
            if isnan(vs(i)) || vs(i) == 0
                qflag_ion_own(i) = -1;
                vs(i) = 0;
            elseif QualityCode_S(i) == 0
                qflag_ion_own(i) = -1;
            end
        end

        % convert for further calculations
        vx_Hz = vx .* 10^6; % [Hz] 
        vs_Hz = vs .* 10^6; % [Hz] 

        % Convert sigma_vx and sigma_vs from MHz to Hz for calculations
        sigma_vx_Hz = sigma_vx .* 10^6; % [Hz] 
        sigma_vs_Hz = sigma_vs .* 10^6; % [Hz] 


        % 3) Ionospheric correction
        % Calculate ionospheric correction for each observation
        iono_x_corr_own = ((vs_Hz.^2) ./ (vx_Hz.^2 - vs_Hz.^2) .* (tau_s - tau_x)); % [s]  
        iono_x_corr_own = iono_x_corr_own';

        % Calculate partial derivatives with respect to each variable
        partial_tau_x = -(vs_Hz.^2) ./ (vs_Hz.^2 - vx_Hz.^2);
        partial_tau_s = (vs_Hz.^2) ./ (vs_Hz.^2 - vx_Hz.^2);
        partial_vx = (tau_s - tau_x) .* (2 * vx_Hz ./ (vs_Hz.^2 - vx_Hz.^2).^2);
        partial_vs = (tau_s - tau_x) .* (2 * vs_Hz .* vx_Hz.^2 ./ (vs_Hz.^2 - vx_Hz.^2).^2);

        % Propagation of uncertainty 
        sigma_iono_x_corr_own = sqrt((partial_tau_x .* sigma_tau_x).^2 + ...
                         (partial_tau_s .* sigma_tau_s).^2 + ...
                         (partial_vx .* sigma_vx_Hz).^2 + ...
                         (partial_vs .* sigma_vs_Hz).^2); % [s]
        sigma_iono_x_corr_own = sigma_iono_x_corr_own';
                     
        % flag where values could not be calculated in the correct way 
        % e.g. if all channel frequencies are missing for one observation
        qflag_ion_own(isnan(iono_x_corr_own) | iono_x_corr_own == 0) = -1;
        qflag_ion_own(isnan(sigma_iono_x_corr_own) | sigma_iono_x_corr_own == 0) = -1;

        % set NaN values to zero, they are already flaged
        iono_x_corr_own(isnan(iono_x_corr_own)) = 0;
        sigma_iono_x_corr_own(isnan(sigma_iono_x_corr_own)) = 0;

        qflag_ion_own = qflag_ion_own';
        sum_qflag_ion_own = sum(qflag_ion_own);

        if sum_qflag_ion_own == 0
            fprintf(' - Same ionospheric delay flag (= 0) used for all scans! \n')
        end

        fprintf('\t Iono corr: directly calculated in VieVS \n')

        % % If you want to filter additional values 
        % msi = median(sigma_iono_x_corr, 'omitnan');
        % siii = 5*msi; % 3 sigma gives too much and 5 sigma too less matching values  

    catch
        fprintf('WARNING vievs_iono: Calculation of the ionospheric correction failed due to a problem of the used data \n')
        iono_x_corr_own = 0;
        sigma_iono_x_corr_own = 0;
        qflag_ion_own = zeros(lenOBS, 1);
        qflag_ion_own(qflag_ion_own == 0) = -1;
    end

end


end