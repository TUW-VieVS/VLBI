% -------------------------------------------------------------------------
%
%                              vex_exper
%
%   Writes the $EXPER Block to a VEX File.
%   General Inforamtion usful for the success of the VLBI administration
%   process.
%
%   Author: 
%       2013-11-12 : Andreas Hellerschmied (heller182@gmx.at)
%   
%   changes       :
%   - 2015-06-16, A. Hellerschmied: Updated for VieVS 2.3
%           
%
%   inputs        :
%   - fid_vex           : File ID of VEX File
%   - vex_para          : Data from VEX paramer file, structure
%   - sched_data        : Scheduling data structure
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
%   - tdays         : Calculate day of year.
%   - invjday
%   
%
%   references    :
%   - VEX File Definition/Example, Rev 1.5b1, 30. Jan. 2002, available
%     online at: http://www.vlbi.org/vex/
%   - VEX Parameter Tables, Rev 1.5b1, 29. Jan. 2002, available
%     online at: http://www.vlbi.org/vex/
%
%-------------------------------------------------------------------------

function [error_code, error_msg] = vex_exper(fid_vex, vex_para, sched_data)

    % Init
    error_code = 0;
    error_msg = '';
    dec_places = 3; % Number of decimal places of the seconds in epoch data.
    
    if (dec_places == 0)
        digits = 2;
    else
        digits = dec_places + 3;
    end
    
    fprintf(fid_vex, '$EXPER;\n');
    fprintf(fid_vex, '*\n');
    
    % Define General Experiment Information:
    try
        fprintf(fid_vex, 'def %s;\n', sched_data.experiment_label);
    catch error
        error_code = 1;
        error_msg = 'Experiment Label not available';  
        return;
    end
    
    try
        fprintf(fid_vex, '    exper_name = %s;\n', sched_data.experiment_label);
    catch error
        error_code = 1;
        error_msg = 'Experiment Label not available';  
        return;
    end
    
    try
        fprintf(fid_vex, '    exper_description = %s;\n', vex_para.exper_description);
    catch error

    end
                
    try
        [year, mon, day, hr, min, sec] = invjday(sched_data.t_nominal_start_jd);
        [days_of_year] = tdays(year, mon, day);
        [year, days_of_year, hr, min, sec] = round_sec_of_date_to_n_dec_figures(year, days_of_year, hr, min, sec, dec_places);
        fprintf(fid_vex, ['    exper_nominal_start = %04.0dy%3sd%02.0fh%02.0fm%0',num2str(digits) ,'.',num2str(dec_places) ,'fs;\n'], year, days_of_year, hr, min, sec);
    catch error
        error_code = 1;
        error_msg = 'Nominal experiment start epoch is not available';  
        return;
    end
    
    try
        [year, mon, day, hr, min, sec] = invjday(sched_data.t_nominal_end_jd);
        [days_of_year] = tdays(year, mon, day);
        [year, days_of_year, hr, min, sec] = round_sec_of_date_to_n_dec_figures(year, days_of_year, hr, min, sec, dec_places);
        fprintf(fid_vex, ['    exper_nominal_stop = %04.0dy%3sd%02.0fh%02.0fm%0',num2str(digits) ,'.',num2str(dec_places) ,'fs;\n'], year, days_of_year, hr, min, sec);
        
        % Account for "delta_t_sec":
%         delta_t_jd = (sched_data.stat(1).obs(sched_data.number_of_obs).delta_t_sec) / (24 * 60 * 60);
%         [year, mon, day, hr, min, sec] = invjday(sched_data.t_nominal_stop_jd + delta_t_jd);
%         [days_of_year] = tdays(year, mon, day);
%         [year, days_of_year, hr, min, sec] = round_sec_of_date_to_n_dec_figures(year, days_of_year, hr, min, sec, dec_places);
%         fprintf(fid_vex, ['    exper_nominal_stop = %04.0dy%3sd%02.0fh%02.0fm%0',num2str(digits) ,'.',num2str(dec_places) ,'fs;\n'], year, days_of_year, hr, min, sec);
    catch error
        error_code = 1;
        error_msg = 'Nominal experiment end epoch is not available';  
        return;
    end
    
    try
        fprintf(fid_vex, '    PI_name = %s;\n', vex_para.PI_name);
    catch error

    end
    
    try
        fprintf(fid_vex, '    PI_email = %s;\n', vex_para.PI_email);
    catch error

    end
    
    try
        fprintf(fid_vex, '    contact_name = %s;\n', vex_para.contact_name);
    catch error

    end
    
    try
        fprintf(fid_vex, '    contact_email = %s;\n', vex_para.contact_email);
    catch error

    end
    
    try
        fprintf(fid_vex, '    scheduler_name = %s;\n', vex_para.scheduler_name);
    catch error

    end
    
    try
        fprintf(fid_vex, '    scheduler_email = %s;\n', vex_para.scheduler_email);
    catch error

    end
    
    try
        fprintf(fid_vex, '    target_correlator = %s;\n', vex_para.target_correlator);
    catch error
        error_code = 1;
        error_msg = 'Target correlator missing.';  
        return;
    end
    
    fprintf(fid_vex, 'enddef;\n');
    
    % Additional notes, if available:
    try
        fprintf(fid_vex, '*    %s\n', vex_para.exper_note);
    catch error

    end

    fprintf(fid_vex, '* ---------------------------------------------------\n');
    fprintf(fid_vex, '*\n');
    
return;

