% -------------------------------------------------------------------------
%
%                              vie_sched_sat_o
%
%   Creates Output for Satellite Scheduling.
%
%   Author: 
%       2013-11-07 : Andreas Hellerschmied
%   
%   changes       :
%   - 2015-06-19, A. Hellerschmied: Updated for VieVS 2.3
%   - 2015-10-19, A. Hellerschmied: Added function call to write combined VEX file 
%   - 2016-05-12, A. Hellerschmied: Bug-fix
%   - 2016-06-14, A. Hellerschmied: - Input argument changed: OUTFILE not needed any more
%                                   - Changes to run this function from "scheduler_interface" (error routines added)
%   - 2016-11-15, A. Hellerschmied: Changes to create combined VEX files in a separate run
%   - 2016-12-22, A. Hellerschmied: Use path to catalog files defind in the VieVS GUI (PARA.pinfile)
%           
%
%   inputs        :
%   - sched_data    : Scheduling data structure
%   - PARA          : Scheduling parameter structure
%     
%
%   outputs       :
%    
%
%   locals        :
% 
%
%   coupling      :
%   - read_parameter_and_cat_files      : Get VEX Parameter from
%                                         VEX-Parameter and -Catalog Files
%   - read_satellite_frequ_file         : Get Satellite Frequency data
%   - get_station_statistics            : Get some staticial data for each observing station
%   - write_vex_file                    : Write VEX file
%   - write_vex_file_comb               : Write combined VEX file
%
%   references    :
%
%-------------------------------------------------------------------------

function [error_code, error_msg] = vie_sched_sat_o(sched_data, PARA)

    % ##### Init #####
    error_code  = 0;
    error_msg   = '';
    out_state   = 1;
    
    % ##### Set file-paths and -names #####
    
    % VEX output dir.:
    pvievs  = '../';
    if PARA.USE_OUTPUT_SUBDIR
        sub_directroy = PARA.OUTPUT_SUBDIR;             % Sub-directory defined in GUI
    else
        sub_directroy = num2str(round(PARA.year_s));    % Sub-directory = Year
    end
    filepath_vex = [pvievs, 'DATA/SCHED/', sub_directroy, '/'];
    if ~isdir(filepath_vex)
        mkdir(filepath_vex);
    end
    
    % Config. and Catalog files:
    filepath_vex_para =     PARA.pinfile;           % From GUI, default: '../CATALOGS/';
    filename_vex_para =     PARA.VEX_PARA_FILENAME; % 'vex_parameter.txt';
    filepath_sat_frequ =    PARA.pinfile;           % From GUI, default: '../CATALOGS/';
    filename_sat_frequ =    PARA.SAT_FREQ_FILENAME; % 'satellite_frequ.txt';
    filepath_vex_cat =      PARA.pinfile;           % From GUI, default: '../CATALOGS/';

        
    while (1)
        
        switch(out_state)

            case 1 % ++++ Get VEX Parameter from VEX-Parameter and -Catalog Files ++++
                [vex_para, error_code, error_msg] = read_parameter_and_cat_files(filepath_vex_para, filename_vex_para, filepath_vex_cat);
                if (error_code == 0)
                    
                    % Save mat-files:
                    save([PARA.pfolder 'vex_para.mat'], 'vex_para');
                    
                    out_state = 2;
                else
                    out_state = 999;
                    error_msg = ['read_parameter_and_cat_files: ', error_msg];
                end

                
            case 2 % ++++ Get Satellite Frequency data ++++
                [sat_frequ, error_code, error_msg] = read_satellite_frequ_file(filepath_sat_frequ, filename_sat_frequ);
                if (error_code == 0)
                    
                    % Save mat-files:
                    save([PARA.pfolder 'sat_frequ.mat'], 'sat_frequ');
                    
                    out_state = 3;
                else
                    out_state = 999;
                    error_msg = ['read_satellite_frequ_file: ', error_msg];
                end
                
                
            case 3 % ++++ Get statistical information ++++
                [sched_data, error_code, error_msg] = get_station_statistics(sched_data);
                if (error_code == 0)
                    % Save "sched_data" mat-files:
                    save([PARA.pfolder 'sched_data.mat'], 'sched_data');
                    switch(PARA.VEX_MODE)
                        case 'combined'
                            out_state = 5;
                        case 'stat_dependent'
                            out_state = 4;
                    end
                else
                    out_state = 999;
                    error_msg = ['get_station_statistics: ', error_msg];
                end
                
                
            case 4 % ++++ Write station dependent VEX file ++++
                [error_code, error_msg] = write_vex_file(sched_data, filepath_vex, vex_para, sat_frequ);
                if (error_code == 0)
                    out_state = 888;
                else
                    out_state = 999;
                    error_msg = ['write_vex_file: ', error_msg];
                end
                
                
            case 5 % ++++ Write combined VEX file ++++
                if ~isfield(vex_para, 'comb_vex_file_ref_station')
                    vex_para.comb_vex_file_ref_station = 'none    ';
                end
                fprintf(' => Reference station defined in vex parameter file: %s\n', vex_para.comb_vex_file_ref_station)
                [error_code, error_msg] = write_vex_file_comb(sched_data, filepath_vex, vex_para, sat_frequ);
                if (error_code == 0)
                    out_state = 888;
                else
                    out_state = 999;
                    error_msg = ['write_vex_file_comb: ', error_msg];
                end

                
            case 888 % ++++ End switch-case routine +++
                fprintf('Creating VEX Files => Finished! \n');
                break; % while loop

                
            case 999 % ++++ ERROR CASE ++++
%                 fprintf('=> ERROR: %s.\n', error_msg);
                break; % while loop
                
        end % switch(out_state)
        
    end % while (1)

return;

