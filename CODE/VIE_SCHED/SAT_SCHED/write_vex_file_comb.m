% -------------------------------------------------------------------------
%
%                              write_vex_file_comb
%
%   This function creatates a combined VEX file for satellite/quasar observations.
%   "Combined" means, that basically the individual station-dependent VEX files are merge, excepting the $SCHED block.
%   This combined VEX file can be used for correlations with DiFX.
%   One station has to be defined for which the information of the $SCHED block is taken.
%
%   Author: 
%       2015-10-19 : Andreas Hellerschmied (andreas.hellerschmied@geo.tuwien.ac.at)
%   
%   changes       :
%   - 2016-06-14, A. Hellerschmied: Input argument changed: filepath instead of filerpath and -name
%   - 2016-11-14, A. Hellerschmied: - Changes to write combined VEX files without stepwise tracking points 
%                                   - Error, if ref. station ID was not found
%           
%
%   inputs        :
%   - sched_data        : Scheduling data structure
%   - filepath_vex      : Filepath of output VEX file
%   - vex_para          : Data from VEX paramer file 
%   - sat_frequ         : Satellite Radio Signal Frequencies
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
%   - vex_header                : Write VEX Heder Lines
%   - vex_global                : Write VEX $GLOBAL Block
%   - vex_exper                 : Write VEX $EXPER Block
%   - vex_station               : Write VEX $STATION Block
%   - vex_procedures            : Write VEX $PROCEDURES Block
%   - vex_site                  : Write VEX $SITE Block
%   - vex_antenna               : Write VEX $ANTENNA Block
%   - vex_das                   : Write VEX $DAS Block
%   - vex_phase_cal_detect      : Write VEX $PHASE_CAL_DETECT Block
%   - vex_if                    : Write VEX $IF Block
%   - vex_bbc                   : Write VEX $BBC Block
%   - vex_tracks                : Write VEX $TRACKS Block
%   - vex_frequ                 : Write VEX $FREQU Block
%   - vex_mode                  : Write VEX $MODE Block
%   - vex_roll                  : Write VEX $ROLL Block
%   - vex_sched                 : Write VEX $SCHED Block
%   - vex_source                : Write VEX $SOURCE Block
%   
%
%   references    :
%   - VEX File Definition/Example, Rev 1.5b1, 30. Jan. 2002, available
%     online at: http://www.vlbi.org/vex/
%   - VEX Parameter Tables, Rev 1.5b1, 29. Jan. 2002, available
%     online at: http://www.vlbi.org/vex/
%
%-------------------------------------------------------------------------

function [error_code, error_msg] = write_vex_file_comb(sched_data, filepath_vex, vex_para, sat_frequ)

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    flag_in_loop = 1; % Main loop flag (if = 0 => exit loop)
    stat_id_list = 1 : length(sched_data.stat);
        
    
    % ##### Prepare reference station used for the combined VEX file #####
    % 8 char station name:
    if isfield(vex_para, 'comb_vex_file_ref_station')
        ref_station_name = sprintf('%-*s' , 8 , vex_para.comb_vex_file_ref_station);
        % Ref. station ID:
        if strcmp(ref_station_name, 'none    ')
           ref_station_id = 0;
           fprintf( ' => Source coordinates (topo.  RaDec) in $SOURCE are taken from first station in scan!\n');
        else
            ref_station_id = find(strncmp({sched_data.stat.name}, ref_station_name, 8));
        end
    else
        ref_station_id = 0;
        fprintf( ' No reference station defined for combined VEX files.\n'); 
        fprintf( ' => Source coordinates (topo.  RaDec) are taken from first station in scan!\n');
    end
        
        

    
    % #### Prepare VEX Filename and Experiment label: ####
    experiment_label = sched_data.exper_name;
    filename_vex = [experiment_label, '_comb.vex'];
    filepathname_new_vex = [filepath_vex, filename_vex];
    sched_data.experiment_label = experiment_label;
    
    
    % ##### loop over all stations #####
    for i_stat = 1 : length(sched_data.stat)
        
        
        % ##### Prepare all required data in "sched_data" structure #####
        
        % Axis type
        if (strcmp(sched_data.stat(i_stat).axis_type, 'AZEL'))
            sched_data.stat(i_stat).axis_type_1 = 'az';
            sched_data.stat(i_stat).axis_type_2 = 'el';
        elseif (strcmp(sched_data.stat(i_stat).axis_type, 'HADC'))
            sched_data.stat(i_stat).axis_type_1 = 'ha';
            sched_data.stat(i_stat).axis_type_2 = 'dec';
        elseif (strcmp(sched_data.stat(i_stat).axis_type, 'XYNS'))
            sched_data.stat(i_stat).axis_type_1 = 'x';
            sched_data.stat(i_stat).axis_type_2 = 'yns'; % N-S - orientation
        elseif (strcmp(sched_data.stat(i_stat).axis_type, 'XYEW'))
            sched_data.stat(i_stat).axis_type_1 = 'x';
            sched_data.stat(i_stat).axis_type_2 = 'yew'; % E-W - orientation
        else
            % ERROR!
            error_code = 1;
            error_msg = 'Unknown or missing axis type.';
            return;
        end
        
        % Def. Labels for VEX $BLOCKS:
        sched_data.stat(i_stat).station_label       = sched_data.stat(i_stat).label;
        sched_data.stat(i_stat).antenna_label       = sched_data.stat(i_stat).name;
        sched_data.stat(i_stat).site_label          = sched_data.stat(i_stat).name;
        
    end % for i_stat = 1 : length(sched_data.stat)
        
        
        

    % ##### open .vex file: ##### 
    fid_vex = fopen(filepathname_new_vex, 'w'); 
    if (fid_vex == -1)
        error_msg = ['Can not open VEX file: ', filepathname_new_vex];
        error_code = 1;
        return;
    end

    
    % #### Wrtite VEX file: ####
    vex_state = 1;

    while(flag_in_loop)

        switch(vex_state)

            
            case 1 % ++++ Header ++++
                % Write note to header:
                vex_para.header_note = ['Combined VEX file for correlator input (Ref. station: ',ref_station_name ,')'];
                [error_code, error_msg] = vex_header(fid_vex, vex_para);
                if (error_code == 0)
                    vex_state = 2;
                else
                    vex_state = 999;
                    error_msg = ['vex_header: ', error_msg];
                end


            case 2 % ++++ $GLOBAL ++++
                [error_code, error_msg] = vex_global(fid_vex, vex_para, sched_data);
                if (error_code == 0)
                    vex_state = 3;
                else
                    vex_state = 999;
                    error_msg = ['vex_global: ', error_msg];
                end


            case 3 % ++++ $EXPER ++++
                [error_code, error_msg] = vex_exper(fid_vex, vex_para, sched_data);
                if (error_code == 0)
                    vex_state = 4;
                else
                    vex_state = 999;
                    error_msg = ['vex_exper: ', error_msg];
                end


            case 4 % ++++ $STATION ++++
                    % [error_code, error_msg] = vex_station(fid_vex, vex_para, sched_data, i_stat);
                    [error_code, error_msg] = vex_station(fid_vex, vex_para, sched_data, stat_id_list);
                    if (error_code == 0)
                        vex_state = 5;
                    else
                        vex_state = 999;
                        error_msg = ['vex_station: ', error_msg];
                    end


            case 5 % ++++ $PROCEDURES ++++
                [error_code, error_msg] = vex_procedures(fid_vex, vex_para, sched_data, stat_id_list);
                if (error_code == 0)
                    vex_state = 6;
                else
                    vex_state = 999;
                    error_msg = ['vex_procedures: ', error_msg];
                end


            case 6 % ++++ $SITE ++++
                [error_code, error_msg] = vex_site(fid_vex, vex_para, sched_data, stat_id_list);
                if (error_code == 0)
                    vex_state = 7;
                else
                    vex_state = 999;
                    error_msg = ['vex_site: ', error_msg];
                end


            case 7 % ++++ $ANTENNA ++++
                [error_code, error_msg] = vex_antenna(fid_vex, vex_para, sched_data, stat_id_list);
                if (error_code == 0)
                    vex_state = 8;
                else
                    vex_state = 999;
                    error_msg = ['vex_antenna: ', error_msg];
                end


            case 8 % ++++ $DAS ++++
                [error_code, error_msg] = vex_das(fid_vex, vex_para, sched_data, stat_id_list);
                if (error_code == 0)
                    vex_state = 9;
                else
                    vex_state = 999;
                    error_msg = ['vex_das: ', error_msg];
                end


            case 9 % ++++ $PHASE_CAL_DETECT ++++
                [error_code, error_msg] = vex_phase_cal_detect(fid_vex, vex_para, sched_data, stat_id_list);
                if (error_code == 0)
                    vex_state = 10;
                else
                    vex_state = 999;
                    error_msg = ['vex_phase_cal_detect: ', error_msg];
                end


            case 10 % ++++ $IF ++++
                [error_code, error_msg] = vex_if(fid_vex, vex_para, sched_data, stat_id_list);
                if (error_code == 0)
                    vex_state = 11;
                else
                    vex_state = 999;
                    error_msg = ['vex_if: ', error_msg];
                end


            case 11 % ++++ $BBC ++++
                [error_code, error_msg] = vex_bbc(fid_vex, vex_para, sched_data, stat_id_list);
                if (error_code == 0)
                    vex_state = 12;
                else
                    vex_state = 999;
                    error_msg = ['vex_bbc: ', error_msg];
                end


            case 12 % ++++ $TRACKS ++++
                [error_code, error_msg] = vex_tracks(fid_vex, vex_para, sched_data, stat_id_list);
                if (error_code == 0)
                    vex_state = 13;
                else
                    vex_state = 999;
                    error_msg = ['vex_tracks: ', error_msg];
                end


            case 13 % ++++ $FREQ ++++
                [error_code, error_msg] = vex_freq(fid_vex, vex_para, sat_frequ, sched_data, stat_id_list);
                if (error_code == 0)
                    vex_state = 14;
                else
                    vex_state = 999;
                    error_msg = ['vex_freq: ', error_msg];
                end


            case 14 % ++++ $MODE ++++
                vex_para.mode_note = ['Ref. station for combined VEX file: ', ref_station_name];
                [error_code, error_msg] = vex_mode(fid_vex, vex_para, sched_data, stat_id_list);
                if (error_code == 0)
                    vex_state = 15;
                else
                    vex_state = 999;
                    error_msg = ['vex_mode: ', error_msg];
                end


            case 15 % ++++ $ROLL ++++
                [error_code, error_msg] = vex_roll(fid_vex, vex_para, sched_data, stat_id_list);
                if (error_code == 0)
                    vex_state = 16;
                else
                    vex_state = 999;
                    error_msg = ['vex_roll: ', error_msg];
                end


            case 16 % ++++ $SCHED ++++
                vex_para.sched_note = ['Ref. station for combined VEX file: ', ref_station_name];
                [error_code, error_msg] = vex_sched(fid_vex, vex_para, sched_data, stat_id_list);
                if (error_code == 0)
                    vex_state = 17;
                else
                    vex_state = 999;
                    error_msg = ['vex_sched: ', error_msg];
                end


            case 17 % ++++ $SOURCE ++++
                % Write a note about the ref. statio name:
                vex_para.source_note = ['Source coordinates are taken from station: ', ref_station_name];
                [error_code, error_msg] = vex_source(fid_vex, vex_para, sched_data, stat_id_list);
                if (error_code == 0)
                    vex_state = 888;
                else
                    vex_state = 999;
                    error_msg = ['vex_source: ', error_msg];
                end


            case 888 % ++++ End switch-case routine ++++
                flag_in_loop = 0;
                fprintf(' => Writing comb. VEX file finished: %s\n', filepathname_new_vex)


            case 999 % ++++ ERROR CASE ++++
                % flag_in_loop = 0;
                % Close pen file:
                if ( (exist('fid_vex', 'var') ) && (fid_vex ~= -1) )
                    fclose(fid_vex);
                end
                return;

        end % switch(vex_state)

    end % while (flag_in_loop)


    
    % ##### close .vex file: #####
    if ( (exist('fid_vex', 'var') ) && (fid_vex ~= -1) )
        fclose(fid_vex);
    end


return;

