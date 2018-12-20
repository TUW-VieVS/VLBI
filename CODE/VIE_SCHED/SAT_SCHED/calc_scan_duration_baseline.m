% #########################################################################
% #     calc_scan_duration_baseline
% #########################################################################
%
% DESCRIPTION
%   This function calculates the duration of a quasar scan for one particular
%   baseline.
%
%
% CREATED  
%   2015-03-31     Andreas Hellerschmied
%
% REFERENCES
% - Jing Sun (2013), VLBI Scheduling strategies with respect to VLBI2010,
%   Geowissenschaftliche Mitteilungen, Heft Nr. 92, ISSN 1811-8380.
%
% COUPLING
% - sobsflux : Calculation of the observed flux density of a source
% - ssefdel
%
%
% INPUT
% - stat_data           - station data structure
% - stat_id_1           - Station IDs of the first station of the treated baseline (referring to the required entry of "stat_data.stat" and "obs_data.stat")
% - stat_id_2           - Station IDs of the second station of the treated baseline (referring to the required entry of "stat_data.stat" and "obs_data.stat")
% - PARA                - Global scheduling parameter strucutre
% - source_quasar       - source (quasars) data structure (ONE sub-set/field of the "source" structure, containing the required data)
% - t_epoch_jd          - Calculation epoch
%
%
% OUTPUT
% - error_code          - Error Code (0 = no erros occured)
% - error_msg           - Error Message (empty, if no errors occured)
% - scan_duration_sec   - Calculated scan duration [sec]
%
% CHANGES:
% - 2015-06-24, A. Hellerschmied: Check for min. scan duration added (PARA.MIN_SCAN)
%

function [scan_duration_sec, error_code, error_msg] = calc_scan_duration_baseline(stat_data, stat_id_1, stat_id_2, el_1, el_2, PARA, source_quasar, t_epoch_jd)
 
    % ##### Init #####
    error_code = 0;
    error_msg = '';

    
    % ##### calculate the observed flux #####

    % Get basline length:
    bl_dx = stat_data.stat(stat_id_1).location.TRF.x - stat_data.stat(stat_id_2).location.TRF.x;
    bl_dy = stat_data.stat(stat_id_1).location.TRF.y - stat_data.stat(stat_id_2).location.TRF.y;
    bl_dz = stat_data.stat(stat_id_1).location.TRF.z - stat_data.stat(stat_id_2).location.TRF.z;
    
    % Get source position:
    ra = source_quasar.ra;
    de = source_quasar.de;
    
    % Get flux parameter:
    fluxpara(1:PARA.MAX_BANDNUM,1:PARA.MAX_FLUXPARA) = source_quasar.fluxpara(1:PARA.MAX_BANDNUM,1:PARA.MAX_FLUXPARA);
    
    % Calc observed flux:
    % Ref.: J. Sun (2013), pp. 8 - 12.
    [obsflux] = sobsflux(t_epoch_jd - 2400000.5, bl_dx, bl_dy, bl_dz, ra, de, fluxpara, PARA);
    
       
    
    % ##### calculate scan duration #####
    
    % Get SEFD parameters:
    sefd_para_1(1:PARA.MAX_BANDNUM,1:PARA.MAX_SEFDPARA) = stat_data.stat(stat_id_1).sefd_para(1:PARA.MAX_BANDNUM,1:PARA.MAX_SEFDPARA);
    sefd_para_2(1:PARA.MAX_BANDNUM,1:PARA.MAX_SEFDPARA) = stat_data.stat(stat_id_2).sefd_para(1:PARA.MAX_BANDNUM,1:PARA.MAX_SEFDPARA);
    
    % Loop init.:
    scan_duration_sec_temp = zeros(PARA.MAX_BANDNUM, 1);
    
    for i_band = 1 : PARA.MAX_BANDNUM % loop over observed bands (1 = X band, 2 = S band)
        
        % #### Calculate adjusted SEFD for both stations by mapping it down to the observed elevation angle ####
        % [sefdel] = ssefdel(el, sefd, y, c0, c1)
        % Ref.: J. Sun (2013), pp. 19 - 20. 
        [sefd_el_1] = ssefdel(el_1, sefd_para_1(i_band,1), sefd_para_1(i_band,2), sefd_para_1(i_band,3), sefd_para_1(i_band,4));
        [sefd_el_2] = ssefdel(el_2, sefd_para_2(i_band,1), sefd_para_2(i_band,2), sefd_para_2(i_band,3), sefd_para_2(i_band,4));
        
        % #### Calculate scan duration for the baseline #### 
        % Ref.: J. Sun (2013), p. 21.
        % Get min. SNR of both stations:
        min_snr = min(stat_data.stat(stat_id_1).min_snr(i_band), stat_data.stat(stat_id_2).min_snr(i_band));
        anum = (1.75 * min_snr / obsflux(i_band)) ^ 2;
        anu1 = sefd_el_1 * sefd_el_2;
        anu2 = PARA.obsmode.samprate * 1.0d6 * PARA.obsmode.chanumband(i_band) * PARA.obsmode.bits;
        scan_duration_sec_temp(i_band) = ceil(anum * (anu1 / anu2) + PARA.CORSYNCH);  
        
    end
    
   
    % ##### Find the max. scan length of the observed bands for the baseline #####
    scan_duration_sec = max(scan_duration_sec_temp(i_band));
    
    % ##### Check min. Scan duration defined in PARA #####
    if scan_duration_sec < PARA.MIN_SCAN
        scan_duration_sec = PARA.MIN_SCAN;
    end
    
    
return;
