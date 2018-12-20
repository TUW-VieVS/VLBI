% #########################################################################
% #     azel_out
% #########################################################################
%
% DESCRITPION
% This function creates AZEL files.
%
% AUTHOR 
%   ???
%
% INPUT
%  - process_list           - Process list with list of sessions
%  - subdir_azel_str        - LEVEL1 sub-driectory to take the data from
%  - path_azel_file_str     - Start time to read in data (optional)
%  - path_azel_file_str     - Output path for Azel files
%  - year                   - year of the session (optional)
%  
% OUTPUT
%
% CHANGES
%  - 2016-08-02, A. Hellerschmied: - Header added
%                                  - Bug fixed (auto. processing did not work after the last update)
%                                  - Handdling of old (VieVS 2.3 or older) and new version (VieVS 3.0) of the sources structure from LEVEL1
%  - 2016-10-13, A. Hellerschmied: Changes applied to support vgosDB and vso files in the process list
%                                   - vgosDB and vso files contain the tags " [vgosDB]" and " [VSO]" respectively in the process list entries. Thsi tag has to be removed to get the actual filename
%
% EXTERNAL CALLS:
% - jul2dat.m
%   


function azel_out(process_list, subdir_azel_str, path_azel_file_str)

% Loop over all input datasets:
for ip = 1 : size(process_list, 1)
    
    yrstr = process_list(ip,1:4);
    sname = process_list(ip,6:end);
    
    % Check, is the session name (in process list) contains the tags " [vgosDB]" or " [VSO]":
    % => This is required to get the actual filename in case of vgosDB or vso files.
    % => Just remove everything after the first blank!
    blank_ind = min(strfind(sname, ' '));
    if ~isempty(blank_ind)
        sname = sname(1 : blank_ind-1);
    end
    
    load(strcat('../DATA/LEVEL1/',subdir_azel_str,'/',num2str(sname),'_scan.mat'));
    load(strcat('../DATA/LEVEL1/',subdir_azel_str,'/',num2str(sname),'_antenna.mat'));
    load(strcat('../DATA/LEVEL1/',subdir_azel_str,'/',num2str(sname),'_sources.mat'));
    
    % create output directory
    if ~exist([path_azel_file_str '/' yrstr], 'dir')
        mkdir([path_azel_file_str '/' yrstr]);
    end
    
    % open output file
    azel_file_path_name = [path_azel_file_str '/' yrstr '/azel_' num2str(sname),'.txt'];
    fid = fopen(azel_file_path_name,'wt');
    
    % write file header
    fprintf(fid,'%%************************************************************\n');
    fprintf(fid,strcat('%% list of az & elevation per observation: ',num2str(sname),'\n'));
    %fprintf(fid,strcat('%% total number of observations: ',num2str(nbas),'\n'));
    fprintf(fid,'%% Columns:\n%%\t 1     .... scannumber\n%%\t 2     .... mjd\n%%\t 3     .... year\n');
    fprintf(fid,'%%\t 4     .... day of year\n%%\t 5-7   .... hour, min, sec\n%%\t 8     .... station\n');
    fprintf(fid,['%%\t 9     .... azimuth [rad]\n%%\t 10    .... elevation' ...
        ' [rad]\n%%\t 11    .... source\n']);
    fprintf(fid,['%%\t 12    .... temperature [deg. C]\n%%\t 13    .... pressure' ...
        ' [hPa]\n']);
    fprintf(fid,'%%\t 14    .... water vapor pressure [hPa]\n');
    fprintf(fid,'%%************************************************************\n');
    fprintf(fid,'%%\n');
    
    % write data to the azel file
    for is = 1 : length(scan)
        mjd = scan(is).mjd;
        
        % #### Get the source name ####
        
        % In case the input source structure has the "new" format (VieVS 3.0) the fiels fields "sources.q" and "sources.s" exist
        % =>  In this case, the field "scan(i_scan).obs_type" is required to distinguish between scans to satellites and to quasars
        if isfield(sources, 'q') || isfield(sources, 's')
            if isfield(scan(is), 'obs_type')
                switch(scan(is).obs_type) 
                    case 'q' % quasar scan
                        source_name_str  = sources.q(scan(is).iso).name;  
                    case 's' % satellite scan
                        source_name_str  = sources.s(scan(is).iso).name; 
                        % Check, if the satellite name has 8 char. 
                        %  => If not => make it 8 char. long! 
                        if length(source_name_str) < 8
                            source_name_str = sprintf('%-*s' , 8 , source_name_str);
                        elseif length(source_name_str) > 8
                            warning('Session %s: Sourcename > 8 char.: %s => %s', sname, source_name_str, source_name_str(1:8));
                            source_name_str = source_name_str(1:8);
                        end
                end
            else
                % Error, if the information on the observation type is missing!
                error('Information on the observation type is missing for session: %s', sname); 
            end
        else
            % "old" format, without the fields "sources.q" and "sources.s"
            source_name_str  = sources(scan(is).iso).name;  
        end
        
        [year,doy,hour,minutes,sec] = jul2dat(mjd);
        for ia = 1:length(scan(is).stat)
            az  = scan(is).stat(ia).az;
            if isempty(az)
            else
                el  = pi/2 - scan(is).stat(ia).zd;
                station = antenna(ia).name; 
                station = sprintf('%-8s', strrep(strtrim(station),' ','_'));   % if there is a blank in between the station name, change it to '_'
                T = scan(is).stat(ia).temp;
                p = scan(is).stat(ia).pres;
                e = scan(is).stat(ia).e;
                
                fprintf(fid, ' %5i %10.5f %4i %3i %2i %2i %5.2f %s %20.15f %20.15f %8s %6.2f %7.2f %6.2f\n',is,mjd,year,doy,hour,minutes,sec,station,az,el,source_name_str,T,p,e);
            end
        end
    end
    
    fclose(fid);
    
    % Status Msg to CW:
    fprintf(1, ' => Created: %s\n', azel_file_path_name);
end

end
