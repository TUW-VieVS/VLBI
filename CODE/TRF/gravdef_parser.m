% Helper function for the creation of the superstationfile with
% mk_superstatFile.m
% Reads the provided gravtiational deformation file, parses it and saves 
% the information in the ns_codes struct. The text file is expected to
% adhear to the following structure:
% 
% # Lines that start with "#" are comments and can be ignored.
% For each antenna the general format is:
%  ANTENNA NUM_PTS SCALE
%  el1  delay1
%  el2  delay2
%  ....
%  ANTEN2 NUM_PTS SCALE
%  EPOCH YYYYMMDD  YYYYMMDD    optional beginning and ending epoch of model
%  el1  delay1
%  el2  delay2
% ...
% Some comments on format.
% 1. ANTENNA is the name of the antenna.
% 2. NUM_PTS are the number of points.
% 3. SCALE is an optional scaling factor. This converts the delay into ps. 
%    If SCALE is 3.336 than the raw values are in mm.   
% 4. Epoch is the optional beginning and ending date for the model. 
% 5. The pairs el1, delay1 are the delay as a function of elevation.
%    The scale factor converts this to ps.
% 
% The delay is saved in ps in the ns_codes.
% 
% Reference for the gravitational deformation: 
%    A complete VLBI delay model for deforming radio telescopes: the Effelsberg case
%    T. Artz · A. Springer · A. Nothnagel
%    J Geod (2014) 88:11451161
%    DOI 10.1007/s00190-014-0749-1
% 

function [ns_codes] = gravdef_parser(ns_codes, gravdefFile)
    fid = fopen(gravdefFile, 'r');
    tline = fgetl(fid);
    
    current_el_delay = [];
    el_delay_idx = 1;
    current_block_info = [];
    epoch_information = [];
    
    while ~feof(fid)
        % try to match a regex for the header line
        header_regex_match = regexp(tline, '^\s*([\w|\-|\_]+)\s+(\d+)\s+(\d+\.\d+)\s*$', 'tokens');
        
        % did it match?
        if ~isempty(header_regex_match)
            % a new header line was found. Is there an old block?
            if ~isempty(current_block_info)
                % find station in ns_codes
                iStat=find(strcmpi(strtrim({ns_codes.name}), current_block_info{1}{1}));

                % if it is empty, try to find stationname with '_' instead of ' '
                if isempty(iStat)
                    iStat=find(strcmpi({ns_codes.name}, ...
                        strrep(current_block_info{1}{1}, ' ', '_')));
                end

                if ~isempty(iStat)
                    if ~isfield(ns_codes(iStat), 'gravdef')
                        curBreak=1;
                    else
                        % get number of break (if there was already the stat found)
                        if isfield(ns_codes(iStat(1)).gravdef, 'break')
                            curBreak=length(ns_codes(iStat).gravdef.break)+1;
                        else
                            curBreak=1;
                        end
                    end
                    iStat=iStat(1);
                    current_el_delay(:,2) = current_el_delay(:,2) * str2double(current_block_info{1}{3});
                    ns_codes(iStat).gravdef.break(curBreak).ez_delay = ...
                        current_el_delay;
                    if ~isempty(epoch_information)
                        ns_codes(iStat).gravdef.break(curBreak).start = ...
                            juliandate(datetime(epoch_information{1}{1},'InputFormat','yyyyMMdd'),'modifiedjuliandate');                       
                        ns_codes(iStat).gravdef.break(curBreak).end   = ...
                            juliandate(datetime(epoch_information{1}{2},'InputFormat','yyyyMMdd'),'modifiedjuliandate');                            
                    else
                        ns_codes(iStat).gravdef.break(curBreak).start = 0;
                        ns_codes(iStat).gravdef.break(curBreak).end   = 99999;
                    end
                else
                    fprintf('%s (station in gravitational deformation file) not found in ns_codes ('' '' -> ''_'' checked)\n', current_block_info{1}{1});
                end
            end
            current_block_info = header_regex_match;
            current_el_delay = zeros(str2double(current_block_info{1}{2}), 2);
            epoch_information = [];
            el_delay_idx = 1;
            
            % advance one line
            tline = fgetl(fid);
            
            % check for optional EPOCH information
            epoch_information = regexp(tline, '^\s*EPOCH\s+(\d{8})\s+(\d{8})\s*$', 'tokens');
            
            
        else
            % try to read el dealy pair
%            el_delay = regexp(tline, '^\s*(\d+\.?\d*)\s+(\-?\d+\.?\d*)\s*$', 'tokens');
            el_delay = sscanf(tline,'%e')';
            if ~isempty(el_delay)
%                current_el_delay(el_delay_idx,:) = str2double(el_delay{1});
                current_el_delay(el_delay_idx,:) = el_delay(1:2);
                el_delay_idx = el_delay_idx + 1;
            end
            
            % advance one line
            tline = fgetl(fid);
        end
    end
    
    % write last block to superstation file
    if ~isempty(current_block_info)
        % find station in ns_codes
        iStat=find(strcmpi(strtrim({ns_codes.name}), current_block_info{1}{1}));

        % if it is empty, try to find stationname with '_' instead of ' '
        if isempty(iStat)
            iStat=find(strcmpi({ns_codes.name}, ...
                strrep(current_block_info{1}{1}, ' ', '_')));
        end

        if ~isempty(iStat)
            if ~isfield(ns_codes(iStat), 'gravdef')
                curBreak=1;
            else
                % get number of break (if there was already the stat found)
                if isfield(ns_codes(iStat(1)).gravdef, 'break')
                    curBreak=length(ns_codes(iStat).gravdef.break)+1;
                else
                    curBreak=1;
                end
            end
            iStat=iStat(1);
            current_el_delay(:,2) = current_el_delay(:,2) * str2double(current_block_info{1}{3});
            ns_codes(iStat).gravdef.break(curBreak).ez_delay = ...
                current_el_delay;
            if ~isempty(epoch_information)
                ns_codes(iStat).gravdef.break(curBreak).start = ...
                    juliandate(datetime(epoch_information{1}{1},'InputFormat','yyyyMMdd'),'modifiedjuliandate');
                ns_codes(iStat).gravdef.break(curBreak).end   = ...
                    juliandate(datetime(epoch_information{1}{2},'InputFormat','yyyyMMdd'),'modifiedjuliandate');
            else
                ns_codes(iStat).gravdef.break(curBreak).start = 0;
                ns_codes(iStat).gravdef.break(curBreak).end   = 99999;
            end
        else
            fprintf('%s (station in gravitational deformation file) not found in ns_codes ('' '' -> ''_'' checked)\n', current_block_info{1}{1});
        end
    end
end