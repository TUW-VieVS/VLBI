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
                            curBreak=length(ns_codes(iStat).ecc.break)+1;
                        else
                            curBreak=1;
                        end
                    end
                    iStat=iStat(1);
                    ns_codes(iStat).gravdef.break(curBreak).ez_delay = ...
                        current_el_delay * str2double(current_block_info{1}{3});
                    if ~isempty(epoch_information)
                        ns_codes(iStat).gravdef.break(curBreak).start = ...
                            mjuliandate(epoch_information{1}{1}, 'yyyymmdd');
                        ns_codes(iStat).gravdef.break(curBreak).end   = ...
                            mjuliandate(epoch_information{1}{2}, 'yyyymmdd');
                    else
                        ns_codes(iStat).gravdef.break(curBreak).start = 0;
                        ns_codes(iStat).gravdef.break(curBreak).end   = 99999;
                    end
                else
                    fprintf('%s (station in ECCDAT.ecc) not found in ns_codes ('' '' -> ''_'' checked)\n', current_block_info{1}{1});
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
            el_delay = regexp(tline, '^\s*(\d+\.?\d*)\s+(\-?\d+\.?\d*)\s*$', 'tokens');
            if ~isempty(el_delay)
                current_el_delay(el_delay_idx,:) = str2double(el_delay{1});
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
                    curBreak=length(ns_codes(iStat).ecc.break)+1;
                else
                    curBreak=1;
                end
            end
            iStat=iStat(1);
            ns_codes(iStat).gravdef.break(curBreak).ez_delay = ...
                current_el_delay * str2double(current_block_info{1}{3});
            if ~isempty(epoch_information)
                ns_codes(iStat).gravdef.break(curBreak).start = ...
                    mjuliandate(epoch_information{1}{1}, 'yyyymmdd');
                ns_codes(iStat).gravdef.break(curBreak).end   = ...
                    mjuliandate(epoch_information{1}{2}, 'yyyymmdd');
            else
                ns_codes(iStat).gravdef.break(curBreak).start = 0;
                ns_codes(iStat).gravdef.break(curBreak).end   = 99999;
            end
        else
            fprintf('%s (station in ECCDAT.ecc) not found in ns_codes ('' '' -> ''_'' checked)\n', current_block_info{1}{1});
        end
    end
end