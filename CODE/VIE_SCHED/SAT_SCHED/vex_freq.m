% -------------------------------------------------------------------------
%
%                              vex_frequ
%
%   Writes the $FREQU Block to a VEX File.
%       - 
%
%   Author: 
%       2013-11-20 : Andreas Hellerschmied (andreas.hellerschmied@geo.tuwien.ac.at)
%   
%   changes       :
%   - 2014-01-13 : Andreas Hellerschmied: 
%       The half bandwith is subtracted from the RF sky frequency to center
%       the satellite carrier frequency within the recorded
%       bandwidth/spectra.
%   - 2015-06-18, A. Hellerschmied: Updated for VieVS 2.3
%   - 2015-06-18, A. Hellerschmied: Added option ("flag_center_sky_frequ") to skip centering the sky frequ. within the recorded band.
%   - 2015-10-20, A. Hellerschmied: Possibility added to write combined VEX
%               files: FREQ blocks are written for all stations defined by the input
%               argument "stat_id_list". Multiple def. labels are skipped.
%   - 2015-11-16, A. Hellerschmied: Minor bug-fix: $FREQ header wwas written once for each station in "stat_id_list".
%   - 2015-11-16, A. Hellerschmied: Option added to use a static freq. setup for satellite scans
%   - 2016-05-12, A. Hellerschmied: Added option to use one static frequ. definition for all satellites ("flag_use_static_freq_setup")
%           
%
%   inputs        :
%   - fid_vex           : File ID of VEX File
%   - vex_para          : Data from VEX paramer file, structure
%   - sat_frequ         : Satellite Frequency structure
%   - sched_data        : scheduling data structure
%   - stat_id_list      : List of IDs of stations which should be considered
%     
%
%   outputs       :
%   - error_code        : Error Code (0 = no erros occured)
%   - error_msg         : Error Message (empty, if no errors occured)
%    
%
%   locals        :
% 
%
%   coupling      :
%   
%
%   references    :
%   - VEX File Definition/Example, Rev 1.5b1, 30. Jan. 2002, available
%     online at: http://www.vlbi.org/vex/
%   - VEX Parameter Tables, Rev 1.5b1, 29. Jan. 2002, available
%     online at: http://www.vlbi.org/vex/
%
%-------------------------------------------------------------------------

function [error_code, error_msg] = vex_freq(fid_vex, vex_para, sat_frequ, sched_data, stat_id_list)

    % Init.:
    error_code = 0;
    error_msg = '';
    temp_str = '';
    flag_first_freq_def = 1;
    all_labels_str_array = {};
    flag_use_static_freq_setup = 0;
    
    % ##### Wite block #####
    if isempty(stat_id_list) % Error
        error_code = 2;
        error_msg = 'Station ID list is empty';  
        return;
    else
        
        % ##### loop over all stations in ID list ##### 
        for i_stat = 1 : length(stat_id_list)
            
            stat_id = stat_id_list(i_stat);
        
            % ##### Get required data #####

            % Station name:
            station_name = sched_data.stat(stat_id).name;

            % Session infos:
            sat_norad_id_vector = sched_data.stat(stat_id).statistics.sat_norad_id_vector;
            sat_name_list       = sched_data.stat(stat_id).statistics.sat_name_list;

            number_of_observed_sat          = length(sched_data.stat(stat_id).statistics.sat_name_list);
            number_of_observed_quasars      = length(sched_data.stat(stat_id).statistics.quasar_name_list);



            % ##### Write $FREQ block to VEX file #####

            % ### Get number of frequ. setup defintions: ###
            % def. 1 => satellite observation setup
            % def. 2 => quasar observation setup
            number_of_freq = eval(['length(vex_para.', station_name , '.freq);']);

            % ### Check, if frequ. setup for quasars is defined, if quasar scans are included: ###
            if number_of_observed_quasars > 0
               if number_of_freq < 2
                    error_code = 1;
                    error_msg = ['Frequency setup for quasar scans is not defined for station: ', station_name];
                    return;
               end
            end
            
            % ### Common frequ. setup for all satellite observations ###
            % ...if flag "flag_one_frequ_def_for_all_sat" is set in vex_freq.cat
            % Use of static RF sky frequ.


            % ### Write the Block headline only once: ###
            if flag_first_freq_def 
                fprintf(fid_vex, '$FREQ;\n');
                fprintf(fid_vex, '*\n');
                flag_first_freq_def = 0;
            end


            % #########################################
            % #### Frequ. def. for satellite scans ####
            % #########################################

            % --- In VEX parameter file (e.g.: vex_parameter.txt) ---
            % ><station_name>.freq(1)= <frequ_def_label_for_satellites>
            
            % Init.:
            i_frequ = 1; % Def. for satellite obs.!
            temp_str = '';
            temp_str_1 = '';
            
            
            % #### Check, if a static frequency setup should be used ####
            % => Only one setup definition in the $FREQ block for all satellites.
            try 
                temp_str_1 = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').flag_use_static_freq_setup.str;']);
                flag_use_static_freq_setup = logical(str2double(temp_str_1));
            catch
                flag_use_static_freq_setup = 0;
            end
            

            if number_of_observed_sat > 0

                % ##### Loop over all observed satellites in this session: #####
                for i_sat = 1 : number_of_observed_sat
                    
                    if flag_use_static_freq_setup && (i_sat == 2)
                        break;
                    end

                    % #### Create satellite label: ####
                    sat_name = sat_name_list{i_sat};
                    sat_label = sat_name(~isspace(sat_name)); % Remove white space


                    % #### Get satelite NORAD ID (as char string): ####
                    sat_number_str = num2str(sat_norad_id_vector(i_sat));
                    
                    
                    % #### Get Sky-Frequencies (frequ_1 and frequ_2) from external source: ####
                    % Init.:
                    sat_frequ_1_MHz = 0;
                    sat_frequ_2_MHz = 0;
                    
                    if ~flag_use_static_freq_setup

                        for i_sat_frequ = 1 : length(sat_frequ.sat)
                            if strncmp(sat_name, sat_frequ.sat(i_sat_frequ).name, 24)

                                % Sky frequency 1:
                                sat_frequ_1_MHz = sat_frequ.sat(i_sat_frequ).frequ_1;

                                % Sky frequency 2, if existing:
                                if isfield(sat_frequ.sat(i_sat_frequ), 'frequ_2')
                                    if ~isempty(sat_frequ.sat(i_sat_frequ).frequ_2)
                                        sat_frequ_2_MHz = sat_frequ.sat(i_sat_frequ).frequ_2;
                                    end
                                end

                                break;
                            end
                        end
                        % Check for availability of satellite carrier frequ:
                        if sat_frequ_1_MHz == 0
                            error_code = 1;
                            error_msg = ['There is no sky frequ. defined for satellite: ', sat_label];  
                            return;
                        end
                    end % if ~flag_use_static_freq_setup


                    % ##### Write to file #####
                    
                    % Init.:
                    temp_str_1 = '';

                    % #### Round sky frequency, if parameter is defined: ####
                    try
                        temp_str_1 = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').round_sky_frequ.str;']);
                        if ~isempty(temp_str_1)
                            round_to_decimals = str2double(deblank(temp_str_1));
                        else
                            round_to_decimals = 0;
                        end
                    catch
                        round_to_decimals = 0;
                    end


                    % #### Center RF sky frequ. in the recorded band: ####
                    try
                        temp_str_1 = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').flag_center_sky_frequ.str;']);
                        if ~isempty(temp_str_1)
                            flag_center_sky_frequ = logical(str2double(deblank(temp_str_1)));
                        else
                            flag_center_sky_frequ = 1;
                        end
                    catch
                        flag_center_sky_frequ = 1;
                    end


                    % --- Definition: ---
                    try
                        temp_str_1 = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').label;']);
                        if ~flag_use_static_freq_setup
                            temp_str = [sat_number_str, '.', temp_str_1];
                        else
                            temp_str = ['sat.', temp_str_1];
                        end


% #####################################
% CHECK LABEL HERE FOR MULTIPLE ENTRIES                
% #####################################
        
                        % Do not write the current def. to the VEX file, if the label was already written to the VEX file!
                        if logical(sum(ismember(all_labels_str_array, temp_str)))
                            continue;
                        end
                        all_labels_str_array = union(all_labels_str_array, temp_str);

                        fprintf(fid_vex, 'def %s;\n', temp_str);
                    catch error
                        error_code = 1;
                        error_msg = 'FREQU Label not available';  
                        return;
                    end


                    % --- Satellite Description: ---
                    try
                        if (sat_frequ_1_MHz ~= 0) && (sat_frequ_2_MHz == 0) % only one sky frequ. defined...
                            fprintf(fid_vex, '* Satellite name: %s, NORAD ID: %s, frequ_1 = %1.4f MHz\n', sat_label, sat_number_str, sat_frequ_1_MHz);
                        elseif (sat_frequ_1_MHz ~= 0) && (sat_frequ_2_MHz ~= 0)% two sky frequ. defined...
                            fprintf(fid_vex, '* Satellite name: %s, NORAD ID: %s, f_1 = %1.4f MHz, f_2 = %1.4f MHz\n', sat_label, sat_number_str, sat_frequ_1_MHz, sat_frequ_2_MHz);
                        elseif ( (sat_frequ_1_MHz == 0) && (sat_frequ_2_MHz == 0) )&& (flag_use_static_freq_setup) % Static freq setup used...
                            fprintf(fid_vex, '* Freq. setup for satellite scans.\n');
                        else
                            error_code = 1;
                            error_msg = 'Error at satellite sky frequ. definition.';  
                            return;
                        end
                    catch error
                        error_code = 1;
                        error_msg = 'Error at satellite sky frequ. definition.';  
                        return;
                    end

                    % --- note, if sky frequency is rounded ---
                    if round_to_decimals ~= 0
                        fprintf(fid_vex, '* RF sky frequencies are rounded to %1.4f MHz\n', round_to_decimals);
                    end


                    % --- freq_note: ---
                    try
                        number_of_notes = eval(['length(vex_para.', station_name , '.freq(', num2str(i_frequ), ').freq_note)']);

                        for i_note = 1 : number_of_notes
                            temp_str = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').freq_note(', num2str(i_note), ').str']);
                            fprintf(fid_vex, '* %s\n', temp_str);
                        end

                    catch error

                    end


                    % --- sample_rate: ---
                    try
                        temp_str = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').sampling_rate.str;']);
                        fprintf(fid_vex, '    sample_rate = %s;\n', temp_str);
                    catch error
                        error_code = 1;
                        error_msg = '"sample_rate" parameter not available';  
                        return;
                    end


                    % --- chan_def: ---
                    number_of_chan_def = eval(['length(vex_para.', station_name , '.freq(', num2str(i_frequ), ').chan_id);']);

                    % Loop over all channels:
                    for i_chan = 1 : number_of_chan_def

                        % Init.:
                        flag_use_sat_sky_frequ_1 = 0;
                        flag_use_sat_sky_frequ_2 = 0;

                        % --- optional parameters ---

                        % band_id:
                        try
                            temp_str = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').band_id(', num2str(i_chan), ').str;']);
                        catch error
                            temp_str = '  ';
                        end


                        % --- necessary parameters: ---


                        % Get channel BW:
                        try
                            str_bbc_chan_bandwidth = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').bbc_chan_bandwidth(', num2str(i_chan), ').str;']);
                            bandwidth_MHz = str2double(str_bbc_chan_bandwidth);
                        catch error
                            error_code = 1;
                            error_msg = '"bbc_chan_bandwidth" parameter not available, or invalid (not a numerical value)';  
                            return;
                        end

                        % rf_sky_freq

                        % ---- Get rf sky frequ.: ----
                        % note: sky frequency in [MHz]
                        % 1.) "rf_sky_freq" does not exist or is empty => ERROR!
                        try
                            temp_str_1 = deblank(eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').rf_sky_freq(', num2str(i_chan), ').str;'])); % Take rf frequ. from the file "vex_freq.cat"
                            if isempty(temp_str_1)
                                error_code = 3;
                                error_msg = ['"rf_sky_freq" parameter not available for station: ', station_name];  
                                return;
                            end
                        catch error
                            error_code = 4;
                            error_msg = ['"rf_sky_freq" parameter not available for station: ', station_name];  
                            %return;
                        end

                        
                        if ~flag_use_static_freq_setup
                            switch(temp_str_1)

                                case 'f1' % 2.) "rf_sky_freq" = "f1" => take frequ_1 value from satellite frequency definition file!
                                    sky_frequ_sat_MHz = sat_frequ_1_MHz;
                                    flag_use_sat_sky_frequ_1 = 1;

                                case 'f2' % 3.) "rf_sky_freq" = "f2" => take frequ_2 value from satellite frequency definition file!
                                    % Check for availability of satellite carrier frequ:
                                    if sat_frequ_2_MHz == 0
                                        error_code = 1;
                                        error_msg = ['There is no sky frequ. 2 defined for the satellite: ', sat_label];  
                                        return;
                                    else
                                        sky_frequ_sat_MHz = sat_frequ_2_MHz;
                                        flag_use_sat_sky_frequ_2 = 1;
                                    end

                                otherwise % 4.) Try if "rf_sky_freq" is available as numerical number [MHz] => take this value!
                                    try
                                        sky_frequ_sat_MHz = str2double(temp_str_1);
                                    catch error
                                        error_code = 5;
                                        error_msg = ['"rf_sky_freq" parameter not available for station: ', station_name];  
                                        return;
                                    end

                            end % switch(temp_str_1)
                        else
                            try
                                sky_frequ_sat_MHz = str2double(temp_str_1);
                            catch error
                                error_code = 5;
                                error_msg = ['"rf_sky_freq" parameter not available for station: ', station_name];  
                                return;
                            end 
                        end % if ~flag_use_static_freq_setup
                        
                        

                        % --- Get frequency offset, if defined in "vex_freq.cat": ---
                        try
                            temp_str_1 = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').rf_sky_freq_offset(', num2str(i_chan), ').str;']); % Take rf frequ. from the file "vex_freq.cat"
                            if ~isempty(temp_str_1)
                                rf_sky_freq_offset_MHz = str2double(deblank(temp_str_1));
                            else
                                rf_sky_freq_offset_MHz = 0;
                            end
                        catch
                            rf_sky_freq_offset_MHz = 0;
                        end

                        % --- Calc. sky frequ. which is written to the VEX file: ---
                        if flag_center_sky_frequ
                            % To "center" the satellite carrier frequency in the middle
                            % of the recorded frequency bandwidth, subtract (bandwith
                            % [MHz)/2 from the RF sky frequency:
                            sky_frequ_out = sky_frequ_sat_MHz - (bandwidth_MHz / 2);
                        else
                            sky_frequ_out = sky_frequ_sat_MHz;
                        end

                        % --- Round sky frequency, if parameter is defined: ---
                        if round_to_decimals ~= 0
                            sky_frequ_out = round(sky_frequ_out / round_to_decimals) * round_to_decimals;
                        end

                        sky_frequ_out = sky_frequ_out + rf_sky_freq_offset_MHz;


        % ############# TEST TEST TEST ################### 
        % Round the sky frequ to 0.1 MHz, necessary for HOBART26:
        % sky_frequ_out = round(sky_frequ_out * 10) / 10;
        % ############# TEST TEST TEST ################### 

                        temp_str_1 = sprintf('%9.4f', sky_frequ_out); 

                        try % Write sky frequ.
                            temp_str = [temp_str, ' : ', temp_str_1, ' MHz'];
                        catch error
                            error_code = 1;
                            error_msg = '"rf_sky_freq" parameter not available';  
                            return;
                        end


                        % sideband_bbc_chan:
                        try
                            temp_str_1 = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').sideband_bbc_chan(', num2str(i_chan), ').str;']);
                            temp_str = [temp_str, ' : ', temp_str_1];
                        catch error
                            error_code = 1;
                            error_msg = '"sideband_bbc_chan" parameter not available';  
                            return;
                        end


                        % write "bbc_chan_bandwidth":
                        temp_str = [temp_str, ' : ', str_bbc_chan_bandwidth, ' MHz'];


                        % chan_id:
                        try
                            temp_str_1 = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').chan_id(', num2str(i_chan), ').str;']);
                            temp_str = [temp_str, ' : ', temp_str_1];
                        catch error
                            error_code = 1;
                            error_msg = '"chan_id" parameter not available';  
                            return;
                        end

                        % bbc_id:
                        try
                            temp_str_1 = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').bbc_id(', num2str(i_chan), ').str;']);
                            temp_str = [temp_str, ' : ', temp_str_1];
                        catch error
                            error_code = 1;
                            error_msg = '"bbc_id" parameter not available';  
                            return;
                        end

                        % phase_cal_id:
                        try
                            temp_str_1 = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').phase_cal_id(', num2str(i_chan), ').str;']);
                            temp_str = [temp_str, ' : ', temp_str_1, ';'];
                        catch error
                            error_code = 1;
                            error_msg = '"phase_cal_id" parameter not available';  
                            return;
                        end


                        % --- Write "chan_def": ---

                        % Create comment:
                        if flag_use_sat_sky_frequ_1 % Frequ_1
                            temp_str_1 = ' * f1';
                        elseif flag_use_sat_sky_frequ_2 % Frequ_2
                            temp_str_1 = ' * f2';
                        else
                            temp_str_1 = '';
                        end

                        % Add the applied RF sky frequency offset as a comment:
                        if rf_sky_freq_offset_MHz ~= 0
                            if isempty(temp_str_1)
                                temp_str_1 = ' * ';
                            else
                                temp_str_1 = [temp_str_1, ', '];
                            end
                            temp_str_1 = [temp_str_1, sprintf('offset = %+1.4f MHz', rf_sky_freq_offset_MHz)];
                        end

                        temp_str = [temp_str, temp_str_1];
                        fprintf(fid_vex, '    chan_def = %s\n', temp_str);



                    end % for i_chan = 1 : number_of_chan_def

                    fprintf(fid_vex, 'enddef;\n');
                    fprintf(fid_vex, '*\n');

                end % for i_sat = 1 : number_of_observed_sat

            end % if number_of_observed_sat > 0



            % ######################################
            % #### Frequ. def. for quasar scans ####
            % ######################################

            % --- In VEX parameter file (e.g.: vex_parameter.txt) ---
            % ><station_name>.freq(2)= <frequ_def_label_for_quasars>


            if number_of_observed_quasars > 0


                % ##### Write to file #####

                 % init:
                temp_str = '';
                temp_str_1 = '';
                i_frequ = 2; % Def. for quasar obs.!


                % #### Round sky frequency, if parameter is defined: ####
                try
                    temp_str_1 = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').round_sky_frequ.str;']);
                    if ~isempty(temp_str_1)
                        round_to_decimals = str2double(deblank(temp_str_1));
                    else
                        round_to_decimals = 0;
                    end
                catch
                    round_to_decimals = 0;
                end


                % #### Center RF sky frequ. in the recorded band: ####
                try
                    temp_str_1 = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').flag_center_sky_frequ.str;']);
                    if ~isempty(temp_str_1)
                        flag_center_sky_frequ = logical(str2double(deblank(temp_str_1)));
                    else
                        flag_center_sky_frequ = 1;
                    end
                catch
                    flag_center_sky_frequ = 1;
                end


                % --- Definition: ---
                try
                    temp_str = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').label;']);


        % #####################################
        % CHECK LABEL HERE FOR MULTIPLE ENTRIES                
        % #####################################
        
                    % Do not write the current def. to the VEX file, if the label was already written to the VEX file!
                    if logical(sum(ismember(all_labels_str_array, temp_str)))
                        continue;
                    end
                    all_labels_str_array = union(all_labels_str_array, temp_str);


                    fprintf(fid_vex, 'def %s;\n', temp_str);
                catch error
                    error_code = 1;
                    error_msg = 'FREQU Label not available';  
                    return;
                end


                % --- Source Description: ---
                fprintf(fid_vex, '* Freu. setup for quasar observations.\n');

                % --- note, if sky frequency is rounded ---
                if round_to_decimals ~= 0
                    fprintf(fid_vex, '* RF sky frequencies are rounded to %1.4f MHz\n', round_to_decimals);
                end

                % --- freq_note: ---
                try
                    number_of_notes = eval(['length(vex_para.', station_name , '.freq(', num2str(i_frequ), ').freq_note)']);
                    for i_note = 1 : number_of_notes
                        temp_str = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').freq_note(', num2str(i_note), ').str']);
                        fprintf(fid_vex, '* %s\n', temp_str);
                    end
                catch error
                end


                % --- sample_rate: ---
                try
                    temp_str = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').sampling_rate.str;']);
                    fprintf(fid_vex, '    sample_rate = %s;\n', temp_str);
                catch error
                    error_code = 1;
                    error_msg = '"sample_rate" parameter not available';  
                    return;
                end


                % --- chan_def: ---
                number_of_chan_def = eval(['length(vex_para.', station_name , '.freq(', num2str(i_frequ), ').chan_id);']);

                % Loop over all channels:
                for i_chan = 1 : number_of_chan_def

                    % --- optional parameters ---

                    % band_id:
                    try
                        temp_str = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').band_id(', num2str(i_chan), ').str;']);
                    catch error
                        temp_str = '  ';
                    end


                    % --- necessary parameters: ---

                    % Get channel BW:
                    try
                        str_bbc_chan_bandwidth = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').bbc_chan_bandwidth(', num2str(i_chan), ').str;']);
                        bandwidth_MHz = str2double(str_bbc_chan_bandwidth);
                    catch error
                        error_code = 1;
                        error_msg = '"bbc_chan_bandwidth" parameter not available, or invalid (not a numerical value)';  
                        return;
                    end

                    % rf_sky_freq

                    % ---- Get rf sky frequ.: ----
                    % note: sky frequency in [MHz]
                    % 1.) "rf_sky_freq" does not exist or is empty => ERROR!
                    try
                        temp_str_1 = deblank(eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').rf_sky_freq(', num2str(i_chan), ').str;'])); % Take rf frequ. from the file "vex_freq.cat"
                        if isempty(temp_str_1)
                            error_code = 3;
                            error_msg = ['"rf_sky_freq" parameter not available for quasar obs. for station: ', station_name];  
                            return;
                        end
                    catch error
                        error_code = 4;
                        error_msg = ['"rf_sky_freq" parameter not available for quasar obs. for station: ', station_name];  
                        %return;
                    end

                    try
                        sky_frequ_quasar_MHz = str2double(temp_str_1);
                    catch error
                        error_code = 5;
                        error_msg = ['"rf_sky_freq" parameter not available for station: ', station_name];  
                        return;
                    end

                    % --- Get frequency offset, if defined in "vex_freq.cat": ---
                    try
                        temp_str_1 = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').rf_sky_freq_offset(', num2str(i_chan), ').str;']); % Take rf frequ. from the file "vex_freq.cat"
                        if ~isempty(temp_str_1)
                            rf_sky_freq_offset_MHz = str2double(deblank(temp_str_1));
                        else
                            rf_sky_freq_offset_MHz = 0;
                        end
                    catch
                        rf_sky_freq_offset_MHz = 0;
                    end


                    % --- Calc. sky frequ. which is written to the VEX file: ---
                    % To "center" the satellite carrier frequency in the middle
                    % of the recorded frequency bandwidth, subtract (bandwith
                    % [MHz)/2 from the RF sky frequency:
                    if flag_center_sky_frequ
                        % To "center" the satellite carrier frequency in the middle
                        % of the recorded frequency bandwidth, subtract (bandwith
                        % [MHz)/2 from the RF sky frequency:
                        sky_frequ_out = sky_frequ_quasar_MHz - (bandwidth_MHz / 2);
                    else
                        sky_frequ_out = sky_frequ_quasar_MHz;
                    end

                    % --- Round sky frequency, if parameter is defined: ---
                    if round_to_decimals ~= 0
                        sky_frequ_out = round(sky_frequ_out / round_to_decimals) * round_to_decimals;
                    end

                    sky_frequ_out = sky_frequ_out + rf_sky_freq_offset_MHz;

                    temp_str_1 = sprintf('%9.4f', sky_frequ_out); 

                    try % Write sky frequ.
                        temp_str = [temp_str, ' : ', temp_str_1, ' MHz'];
                    catch error
                        error_code = 1;
                        error_msg = '"rf_sky_freq" parameter not available';  
                        return;
                    end


                    % sideband_bbc_chan:
                    try
                        temp_str_1 = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').sideband_bbc_chan(', num2str(i_chan), ').str;']);
                        temp_str = [temp_str, ' : ', temp_str_1];
                    catch error
                        error_code = 1;
                        error_msg = '"sideband_bbc_chan" parameter not available';  
                        return;
                    end

                    % write "bbc_chan_bandwidth":
                    temp_str = [temp_str, ' : ', str_bbc_chan_bandwidth, ' MHz'];


                    % chan_id:
                    try
                        temp_str_1 = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').chan_id(', num2str(i_chan), ').str;']);
                        temp_str = [temp_str, ' : ', temp_str_1];
                    catch error
                        error_code = 1;
                        error_msg = '"chan_id" parameter not available';  
                        return;
                    end

                    % bbc_id:
                    try
                        temp_str_1 = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').bbc_id(', num2str(i_chan), ').str;']);
                        temp_str = [temp_str, ' : ', temp_str_1];
                    catch error
                        error_code = 1;
                        error_msg = '"bbc_id" parameter not available';  
                        return;
                    end

                    % phase_cal_id:
                    try
                        temp_str_1 = eval(['vex_para.', station_name , '.freq(', num2str(i_frequ), ').phase_cal_id(', num2str(i_chan), ').str;']);
                        temp_str = [temp_str, ' : ', temp_str_1];
                    catch error
                        error_code = 1;
                        error_msg = '"phase_cal_id" parameter not available';  
                        return;
                    end


                    % --- Write "chan_def": ---

                    % Add the applied RF sky frequency offset as a comment:
                    if rf_sky_freq_offset_MHz ~= 0

                        temp_str_1 = sprintf('; * offset = %+1.4f MHz', rf_sky_freq_offset_MHz);
                    else
                        temp_str_1 = ';';
                    end

                    temp_str = [temp_str, temp_str_1];
                    fprintf(fid_vex, '    chan_def = %s\n', temp_str);


                end % for i_chan = 1 : number_of_chan_def

                fprintf(fid_vex, 'enddef;\n');
                fprintf(fid_vex, '*\n');

            end % if number_of_observed_quasars > 0
            
        end % for i_stat = 1 : length(stat_id_list)
    
        % #### Additional note, if available: ####
        try
            fprintf(fid_vex, '*    %s\n', vex_para.freq_note);
        catch error
        end

        fprintf(fid_vex, '* ---------------------------------------------------\n');
        fprintf(fid_vex, '*\n');
        
    end % if isempty(stat_id_list)
        
return;
    
    
