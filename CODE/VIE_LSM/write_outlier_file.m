% #########################################################################
% #     write_outlier_file
% #########################################################################
%
% DESCRITPION
% This function writes the outliers detected in VIE_LSM to an ASCII file.
%
% Format:
% | Station name 1 (8 char) | blank | Station name 2 (8 char) | blank | MJD
%
% NOTE:
% The following fields are requird in the parameter structure (general session info):
%  - parameter.session_name
%  - parameter.year
% => This might cause problems with LEVEL 1 data from older VieVS versions (before V3.0)!
%
%
% AUTHOR
%       A. Hellerschmied (2016-10-12)
%
% INPUT
%   out_v           Vector with outlier indics (outliers in "oc_observ")
%   antenna         VieVS antenna structure
%   scan            VieVS scan structure
%   parameter       VieVS parameter structure
%   
%
% OUTPUT
%   out             Matrix containing the scan IDs (out(i,1)) and observation IDs (out(i,2)) of outliers
%
% CHANGES
%

function out = write_outlier_file(out_v, antenna, scan, parameter, sources)

if ~isempty(out_v)
    % ##### Init.: #####
    out     = [];
    inum    = 0;
    ico     = 0;
    
    % ##### Get sattion names and MJD epochs for each outlier in "out_v" #####
    for isc = 1 : length(scan)
        for iobs = 1 : length(scan(isc).obs)
            inum = inum + 1;
            for iout = 1 : size(out_v,2)
                if inum == out_v(iout)
                    ico = ico + 1;
                    stat1_out(ico).name = antenna(scan(isc).obs(iobs).i1).name;
                    stat2_out(ico).name = antenna(scan(isc).obs(iobs).i2).name;
                    mjd_out(ico) = scan(isc).mjd;
                    out(ico,1) = isc; 
                    out(ico,2) = iobs;
                    sou(ico).name = sources.q(scan(isc).iso).name;
                end
            end
        end
    end
    
    % ##### File directory and path #####
    
    % Get output dir.:
    if isfield(parameter, 'year')
        if ~isempty([parameter.lsmopt.dirout])
            outdir = ['../DATA/OUTLIER/', parameter.lsmopt.dirout, '/', parameter.year, '/'];
        else
            outdir = ['../DATA/OUTLIER/', parameter.year, '/'];
        end
    else
        % If "parameter.year" is not availabe, the outlier file is written to /DATA/OUTLIER/
        outdir = '../DATA/OUTLIER/';
    end
    
    % Check output dir.:
    if ~exist(outdir, 'dir')
        mkdir(outdir);
        fprintf('Creating new directory for outlier files: %s\n',outdir);
    end
    
    % Get filename:
    if isfield(parameter, 'session_name')
        session_name = parameter.session_name;  % Since VieVS 3.0
    else
        session_name = antenna(1).session;      % Older VieVS versions
    end
   
    % Full filepath:
    out_name_path_str = [outdir, session_name, '.OUT'];
    
    % Status msg.:
    if exist(out_name_path_str, 'file')
        fprintf('Outliers are appended to existing outlier file: %s\n', out_name_path_str);
    else
        fprintf('Outliers are written to new outlier file: %s\n', out_name_path_str);
    end
       
    % Write to file:
    fid_out = fopen(out_name_path_str,'a');
    if fid_out == -1
        error('Outlier file (%s) could not be openend!', out_name_path_str);
    end
    for isc_out = 1 : size(out, 1)
        fprintf(fid_out, '%s %s %4.12f %s\n' ,stat1_out(isc_out).name, stat2_out(isc_out).name, mjd_out(isc_out), sou(isc_out).name);
    end
    fclose(fid_out);
    
    % Notification:
    fprintf('The %d detected outlier observations are NOT eleminated!\n', length(out_v));
end

