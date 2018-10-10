% #########################################################################
% #     check_scan_duration_and_calc_repos_epochs
% #########################################################################
%
% DESCRIPTION
%   This function checks, if the entered scan period (from t_start_jd to t_end_jd) 
%   is equal (up to a certain treshold defined by: PARA.DELTA_REPOS_INT) to a 
%   integer number of reposition intervals (defined by PARA.REPOS_INT).
%   
%   It the entered duration is OK:
%    - Calculate the antenna repositioning epochs (repos_epochs_jd)
%
%   It the entered duration is NOT OK:
%    - Calculate an alternative scan end time (alt_t_end_jd)
%   
%   Required for stepwise satellite tracking!
%
%
%
% CREATED  
%   2015-04-02     Andreas Hellerschmied
%
% REFERENCES
%
% COUPLING
%
%
% INPUT
% - PARA                - Global scheduling parameter strucutre
% - t_start_jd          - Scan start time [JD]
% - t_end_jd            - Scan end time [JD]
%
% OUTPUT
% - error_code              - Error Code (0 = no erros occured)
% - error_msg               - Error Message (empty, if no errors occured)
% - flag_scan_duration_ok   - 1, if the scan duration is
% - alt_t_end_jd            - Alternative scan end time [JD], if the entered one is not OK
% - repos_epochs_jd         - Antenna reposition epochs, calculated for the current case (interval: PARA.REPOS_INT, thrashold: PARA.DELTA_REPOS_INT)
%
% CHANGES:
%   - 2016-08-02: A. Hellerschmied: Ther was a bug when checking the min. scan duration. This is fixed now.
%

function [flag_scan_duration_ok, alt_t_end_jd, repos_epochs_jd, error_code, error_msg] = check_scan_duration_and_calc_repos_epochs(PARA, t_start_jd, t_end_jd)
 
    % ##### Init #####
    error_code = 0;
    error_msg = '';
    flag_scan_duration_ok = 0;
    alt_t_end_jd = [];
    repos_epochs_jd = [];
    

    % ##### Check scan duration #####
    
    % Calc. scan duration
    scan_duration_sec_temp = (t_end_jd - t_start_jd) * 86400;
    
    % Round to milli-seconds, to get rid of numerical problems...
    scan_duration_sec_temp = round(scan_duration_sec_temp*1000)/1000;
    
    % Check mimimum scan length:
    if scan_duration_sec_temp < PARA.MIN_SAT_SCAN_TIME
        fprintf(' => Scan duration (%3.1f) is smaller than the min. scan duration (PARA.MIN_SAT_SCAN_TIME = %3.1f).\n', scan_duration_sec_temp, PARA.MIN_SAT_SCAN_TIME);
        flag_scan_duration_ok = 0;
    else
        if mod(scan_duration_sec_temp, PARA.REPOS_INT) > PARA.DELTA_REPOS_INT % the scan duration is not a multiple of the antenna repos. interval...
            flag_scan_duration_ok = 0;
        else % Scan duration is OK...
            flag_scan_duration_ok = 1;
        end
    end

    % ##### Entered scan duration is NOT OK #####
    
    % Calc. alternative scan end time:
    if ~flag_scan_duration_ok
        multiples = floor(scan_duration_sec_temp / PARA.REPOS_INT);
        if multiples == 0 % If the scan duration is shorter than "PARA.REPOS_INT"...
            multiples = 1;
        end
        scan_duration_sec_temp = multiples * PARA.REPOS_INT;
        alt_t_end_jd = t_start_jd + scan_duration_sec_temp / 86400;
        
        
    % ##### Entered scan duration is OK #####    
    elseif flag_scan_duration_ok
    
        % Number of antenna repos epochs:
%         number_of_epochs = round(scan_duration_sec_temp / PARA.REPOS_INT) + 1;
        number_of_epochs = ceil(scan_duration_sec_temp / PARA.REPOS_INT) + 1;
    
        % Calc. antenna repos. epochs:
        
        t_epoch_jd_temp = t_start_jd;
        repos_epochs_jd = zeros(number_of_epochs, 1);
        
        for i_epoch = 1 : (number_of_epochs - 1)
            repos_epochs_jd(i_epoch) = t_epoch_jd_temp;
            t_epoch_jd_temp = t_epoch_jd_temp + (PARA.REPOS_INT / (24*60*60) ); % (24*60*60) => Conversion from [sec] to [days]
        end
        
        % Assign session end time as end of last repos. interval:
        repos_epochs_jd(number_of_epochs) = t_end_jd;
        
        
    end
    

return;
