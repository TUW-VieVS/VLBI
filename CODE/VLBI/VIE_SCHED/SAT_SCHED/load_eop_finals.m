% #########################################################################
% #     load_eop_finals
% #########################################################################
%
% DESCRITPION
% This function loads EOP finals data from VieVS specific EOP data files for the
% defined timer period.
%
% AUTHOR 
%   A. Hellerschmied, 18.02.2015
%
% INPUT
%   path_eop        EOP data directory
%   filename_eop    EOP file name
%   mjd_start       Begion of period to get EOP data for [MJD]
%   mjd_end         End of period to get EOP data for [MJD]
%
% OUTPUT
%   mjd             MJD of EOP data epoch [d]
%   xp_rad          Polar motion [rad]
%   yp_rad          Polar motion [rad]
%   dut1_sec        UT1 - UTC [sec]
%   dX_rad          Celestial pole offset [rad]
%   dY_rad          Celestial pole offset [rad]
%   error_code  Error Code (0 = no erros occured)
%   error_msg   Error Message (empty, if no errors occured)
%
% CHANGES
% - 2017-03-21, A. Hellerschmied: Now the EOP finals file directly downloaded from the IERS webpage can be loaded.
%   

function [mjd, xp_rad, yp_rad, dut1_sec, dX_rad, dY_rad, error_code, error_msg] = load_eop_finals(path_eop, filename_eop, mjd_start, mjd_end)

    % init.:
    error_code = 0;
    error_msg = '';
    mjd = 0;
    xp = 0;
    yp = 0;
    dut1 = 0;
    dX = 0;
    dY = 0;
    
    period_extension = 1; % Extension of input period for providing EOP data [d]
    
    % const.:
    mas2rad_const = (1/1000) * (pi/180) * (1/3600);
    
    % ##### Get EOP epochs #####
    mjd_min = floor(mjd_start) - period_extension;
    mjd_max = ceil(mjd_end) + period_extension;
    
    
    
    % ##### Load data from EOP files #####
    % Open file:
    eop_file_str = [path_eop, filename_eop];
    fid = fopen(eop_file_str);
    if fid == -1
        error_code = 1;
        error_msg = ['Failed to open EOP file: ', eop_file_str];  
        return;
    else
        fprintf(1, '   - Loading EOP file: %s\n', eop_file_str);
    end

    % Read data:
    eop_data = textscan(fid, '%2c%2c%2c%1c%8c%1c%1c%1c%9c%9c%1c%9c%9c%2c%1c%10c%10c%1c%7c%7c%2c%1c%1c%9c%9c%1c%9c%9c%10c%10c%11c%10c%10c%*[^\n]', 'Whitespace', '');

    % Close file:
    fclose(fid);

    % Get epoch indices:
    mjd = str2double(cellstr(eop_data{5}));
    mjd_ind = (mjd >= mjd_min) & (mjd <= mjd_max);

    % Check, if EOP data is available for the observation epochs:
    if sum(mjd >= mjd_max) == 0
        error_code = 2;
        error_msg = ['No EOP data available for scheduled epochs (+-', sprintf('%d', period_extension), ' days for interpolation)! EOP file: ', eop_file_str];  
        return;
    end

    % Prepare data vectors:
    mjd = mjd(mjd_ind);
    tmp_cellstr = cellstr(eop_data{9});
    xp_arcsec = str2double(tmp_cellstr(mjd_ind));
    xp_arcsec(isnan(xp_arcsec)) = 0;
    xp_mas = xp_arcsec * 1000;

    tmp_cellstr = cellstr(eop_data{12});
    yp_arcsec = str2double(tmp_cellstr);
    yp_arcsec = yp_arcsec(mjd_ind);
    yp_arcsec(isnan(yp_arcsec)) = 0;
    yp_mas = yp_arcsec * 1000;

    tmp_cellstr = cellstr(eop_data{16});
    dut1_sec = str2double(tmp_cellstr);
    dut1_sec = dut1_sec(mjd_ind); 
    dut1_sec(isnan(dut1_sec)) = 0;
    dut1_msec = dut1_sec * 1000;

    tmp_cellstr = cellstr(eop_data{24});
    dX_mas = str2double(tmp_cellstr);
    dX_mas = dX_mas(mjd_ind);
    dX_mas(isnan(dX_mas)) = 0;

    tmp_cellstr = cellstr(eop_data{27});
    dY_mas = str2double(tmp_cellstr);
    dY_mas = dY_mas(mjd_ind); 
    dY_mas(isnan(dY_mas)) = 0;

    % Unit conversions:
    xp_rad = xp_mas * mas2rad_const;    % [rad]
    yp_rad = yp_mas * mas2rad_const;    % [rad]
    dX_rad = dX_mas * mas2rad_const;    % [rad]
    dY_rad = dY_mas * mas2rad_const;    % [rad]
    dut1_sec = dut1_msec / 1000;        % [sec]

return;