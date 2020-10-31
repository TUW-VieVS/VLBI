% ************************************************************************
%   Description:
%   This function reads the OPT files. 
%   
%   Additional Options: You can specify a list of stations to be removed
%   from every session in your list. This is done in the code. Set
%   remove_sprecial_stations to true and write your station list in the
%   cell array stations_to_be_removed
%
%   Input:	
%      OPT filename.
%
%   Output:
%   - ini_opt
%   - bas_excl
% 
%   External calls: 	
%       
%   Coded for VieVS: 
%   Jul 2012 by Matthias Madzak
%
%   Revision: 
%   11 Jan 2014 by Lucia Plank: down time for excluded station added (11 Feb 2014, corr. by Hana) 
% 	03 Apr 2014 by Matthias Madzak: Added clock breaks to be read by this function
%   16 May 2014 by Monika Tercjak adding reading new field "Stations to be down-weighted" 
%   05 Oct 2015 by A. Hellerschmied: - Parsing of parameters changed: From now on it doesn't matter, how many blanks are included within one line. Only the order of the arguments is relevant.
%                                    - Added possibility to enter MJD values (start, stop) to exclude data for specified time spans ("STATIONS TO BE EXCLUDED")
%   04 Dec 2015 by A. Hellerschmied: - Clock breaks are also loaded here now.
%                                    - Added possibility to enter clock break epochs in the following format: YYMMDDhhmmss
%   21 Dec 2015 by A. Hellerschmied: Remove clock breaks of excluded stations from ini_opt.clk_break struct (...also clk. breaks during excluded time periods)
%   07 Jan 2016 by A. Hellerschmied: - Matlab backward compatibility problem solved (textscan, multiple delimiters)
%                                    - Do not remove clock breaks at stations, if only a subset of the observation data is excluded.
%   28 Jan 2016 by C. Sch�nberger: '-' are allowed in station name. sources can be excluded for a certain time span.
%   02 Feb 2016 by A. Girdiuk: bug-fix in log-messages
%   09 Feb 2016 by C. Sch�nberger: Bug-fix: Date/time input (YYDDMMhhmm) now works for 19xx and 20xx 
%   11 Feb 2016 by D. Mayer: It is now possible to append comments (beginning with #) to OPT statements
%   31 Mar 2017 by D. Mayer: added the possibility to remove list of station from every session in the code
%   28 Aug 2018 by D. Landskron: bug correct corrected with excluding only a time frame of a station
% ************************************************************************

function [ini_opt, bas_excl]=readOPT(optfil,remove_sprecial_stations,stations_to_be_removed)

station_remove_flags = true(size(stations_to_be_removed,1),1);

% #### Options ####
blank_replace_char = '_'; % blanks in station names are replaced by this character.

% ##### Init. #####
ini_opt.sta_excl='';
ini_opt.refclock='';
ini_opt.sour_excl='';
ini_opt.stat_dw=''; % Monika
ini_opt.bas_excl=[];
ini_opt.scan_excl=[];
ini_opt.no_cab='';
ini_opt.num_clk_breaks = 0;
ini_opt.clk_break.stat_name = [];
ini_opt.clk_break.mjd = [];
ini_opt.bdco_est.sta1='';
ini_opt.bdco_est.sta2='';

bas_excl='';


% ##### Open file #####
fid = fopen(optfil,'r');
  

% ##### Loop over all lines #####
while ~feof(fid)
	
    % #### Read one line ####
    str = fgetl(fid);
    
    
    % ##### Parse parameter blocks: #####

    % ### STATIONS TO BE EXCLUDED ###
    if strcmp(str(1:min(length(str),24)),'STATIONS TO BE EXCLUDED:')
        temp_str = textscan(str, '%s', 'CommentStyle', '#'); % parse line
        multiples = str2double(temp_str{1}{end});
        for nex = 1 : multiples
            str = fgetl(fid); % Get line
            
            flag_excl_station = false;
            flag_excl_time_period = false;
            temp_str = textscan(str, '%s', 'delimiter', ' -', 'MultipleDelimsAsOne',1, 'CommentStyle', '#'); % parse line, delimiters: ' ', '-'
            
            if remove_sprecial_stations
                index_special_removed_station = strcmp(temp_str{1}{1},stations_to_be_removed);
                if sum(index_special_removed_station) > 0 
                    temp_str{1} = stations_to_be_removed(index_special_removed_station);
                    station_remove_flags(index_special_removed_station) = 0;
                end
            end
            
            switch size(temp_str{1}, 1)
                case 1 % Exlude station
                    station_name_str = temp_str{1}{1};
                    flag_excl_station = true; 
                case 2 % Exlude station, 8 char. station name with blank or -
                    [flag_ok, station_name_str] = checkStatNameArrayFields(temp_str, 1, 2, str, blank_replace_char);
                    flag_excl_station = true;
                case 3 % % Exlude time period
                    station_name_str = temp_str{1}{1};
                    start_str = temp_str{1}{2};
                    end_str = temp_str{1}{3};
                    flag_excl_time_period = true;
                case 4 % Exlude time period, 8 char. station name with blank or -
                    [flag_ok, station_name_str] = checkStatNameArrayFields(temp_str, 1, 2, str, blank_replace_char);
                    start_str = temp_str{1}{3};
                    end_str = temp_str{1}{4};
                    flag_excl_time_period = true;
                otherwise  % Error case
                    error('STATIONS TO BE EXCLUDED: Invalid input argument.');			%%% 2016-02-02
            end
            % Check station name length:
            [flag_ok, station_name_str] = checkStatNameLength(station_name_str);
            ini_opt.sta_excl(nex,:) = station_name_str;
            if flag_excl_time_period
                if (length(start_str) == 10) && isempty(strfind(start_str, '.')) % Format: YYMMDDhhmm
                    if str2double(start_str(1:2))<79
                        ini_opt.sta_excl_start(nex) = modjuldat(str2double(['20', start_str(1:2)]), str2double(start_str(3:4)), str2double(start_str(5:6)), str2double(start_str(7:8)), str2double(start_str(9:10)));
                    else
                        ini_opt.sta_excl_start(nex) = modjuldat(str2double(['19', start_str(1:2)]), str2double(start_str(3:4)), str2double(start_str(5:6)), str2double(start_str(7:8)), str2double(start_str(9:10)));
                    end
                else % MJD
                    ini_opt.sta_excl_start(nex) = str2double(start_str);
                end
                if (length(end_str) == 10) && isempty(strfind(end_str, '.')) % Format: YYMMDDhhmm
                    if str2double(end_str(1:2))<79
                        ini_opt.sta_excl_end(nex) = modjuldat(str2double(['20', end_str(1:2)]), str2double(end_str(3:4)), str2double(end_str(5:6)), str2double(end_str(7:8)), str2double(end_str(9:10)));
                    else
                        ini_opt.sta_excl_end(nex) = modjuldat(str2double(['19', end_str(1:2)]), str2double(end_str(3:4)), str2double(end_str(5:6)), str2double(end_str(7:8)), str2double(end_str(9:10)));
                    end
                else % MJD
                    ini_opt.sta_excl_end(nex) = str2double(end_str);
                end
            elseif flag_excl_station
                ini_opt.sta_excl_start(nex) = 0;
                ini_opt.sta_excl_end(nex) = 0;
            else
                error('Exclude stations: Unknown mode.');
            end
        end
        
        
    % ### CLOCK BREAKS ###    
    elseif strcmpi(str(1:min(length(str),13)), 'CLOCK BREAKS:')
        temp_str = textscan(str, '%s', 'CommentStyle', '#'); % parse line
        ini_opt.num_clk_breaks = str2double(temp_str{1}{end});
        for nex = 1 : ini_opt.num_clk_breaks
            str = fgetl(fid); % Get line
            temp_str = textscan(str, '%s', 'CommentStyle', '#'); % parse line
            switch size(temp_str{1}, 1)
                case 2
                    station_name_str = temp_str{1}{1};
                    epoch_str = temp_str{1}{2};
                case 3
                    [flag_ok, station_name_str] = checkStatNameArrayFields(temp_str, 1, 2, str, blank_replace_char);
                    epoch_str = temp_str{1}{3};
                otherwise  % Error case
                    error('CLOCK BREAKS: Invalid input argument.');					%%% 2016-02-02
            end
            [flag_ok, station_name_str] = checkStatNameLength(station_name_str);
            ini_opt.clk_break(nex,:).stat_name = station_name_str;
            if (length(epoch_str) == 12) && isempty(strfind(epoch_str, '.')) % Format: YYMMDDhhmmss
                if str2double(epoch_str(1:2))<79
                    ini_opt.clk_break(nex,:).mjd = modjuldat(str2double(['20', epoch_str(1:2)]), str2double(epoch_str(3:4)), str2double(epoch_str(5:6)), str2double(epoch_str(7:8)), str2double(epoch_str(9:10)), str2double(epoch_str(11:12)));
                else
                    ini_opt.clk_break(nex,:).mjd = modjuldat(str2double(['19', epoch_str(1:2)]), str2double(epoch_str(3:4)), str2double(epoch_str(5:6)), str2double(epoch_str(7:8)), str2double(epoch_str(9:10)), str2double(epoch_str(11:12)));
                end
            else % MJD
                ini_opt.clk_break(nex,:).mjd = str2double(epoch_str);
            end
        end
        
        
    % ### CLOCK REFERENCE ###
    elseif strcmpi(str(1:min(length(str),16)), 'CLOCK REFERENCE:')
        str = fgetl(fid); % read line
        temp_str = textscan(str, '%s', 'CommentStyle', '#'); % parse line
        switch size(temp_str{1}, 1)
            case 1
                station_name_str = temp_str{1}{1};
            case 2
                [flag_ok, station_name_str] = checkStatNameArrayFields(temp_str, 1, 2, str, blank_replace_char);
            otherwise  % Error case
                error('CLOCK REFERENCE: Invalid input argument.');
        end
        [flag_ok, station_name_str] = checkStatNameLength(station_name_str);
        ini_opt.refclock(1,:) = station_name_str;


    % ### SOURCES TO BE EXCLUDED ###
    elseif strcmp (str(1:min(length(str),23)),'SOURCES TO BE EXCLUDED:')
        temp_str = textscan(str, '%s', 'CommentStyle', '#'); % parse line
        multiples = str2double(temp_str{1}{end});
        for nex = 1 : multiples
            minus=0;
            str = fgetl(fid); % Get line
            temp_str = textscan(str, '%s', 'delimiter', ' -', 'MultipleDelimsAsOne', 1, 'CommentStyle', '#');% parse line, delimiters: ' ', '-'
            source_name_str = temp_str{1}{1};
            if mod(length(temp_str{1,1}),2)==0
                source_name_str=strcat(temp_str{1}{1},'-',temp_str{1}{2});
                minus=1;
            end
            source_name_str = sprintf('%-*s' , 8 , source_name_str); % 8 char source name => add blanks, if necessary 
            ini_opt.sour_excl(nex,:) = source_name_str;
            if length(temp_str{1,1})>2
                if minus==1
                    start_str = temp_str{1}{3};
                    end_str = temp_str{1}{4};
                else
                    start_str = temp_str{1}{2};
                    end_str = temp_str{1}{3};
                end
                if (length(start_str) == 10) && isempty(strfind(start_str, '.')) % Format: YYMMDDhhmm
                    if str2double(start_str(1:2))<79
                        ini_opt.sour_excl_start(nex) = modjuldat(str2double(['20', start_str(1:2)]), str2double(start_str(3:4)), str2double(start_str(5:6)), str2double(start_str(7:8)), str2double(start_str(9:10)));
                    else
                        ini_opt.sour_excl_start(nex) = modjuldat(str2double(['19', start_str(1:2)]), str2double(start_str(3:4)), str2double(start_str(5:6)), str2double(start_str(7:8)), str2double(start_str(9:10)));
                    end
                else % MJD
                    ini_opt.sour_excl_start(nex) = str2double(start_str);
                end
                if (length(end_str) == 10) && isempty(strfind(end_str, '.')) % Format: YYMMDDhhmm
                    if str2double(end_str(1:2))<79
                        ini_opt.sour_excl_end(nex) = modjuldat(str2double(['20', end_str(1:2)]), str2double(end_str(3:4)), str2double(end_str(5:6)), str2double(end_str(7:8)), str2double(end_str(9:10)));
                    else
                        ini_opt.sour_excl_end(nex) = modjuldat(str2double(['19', end_str(1:2)]), str2double(end_str(3:4)), str2double(end_str(5:6)), str2double(end_str(7:8)), str2double(end_str(9:10)));
                    end
                else % MJD
                    ini_opt.sour_excl_end(nex) = str2double(end_str);
                end
             else % Exclude the complete time span
                ini_opt.sour_excl_start(nex) = 0;
                ini_opt.sour_excl_end(nex) = 0;
            end                              
        end

    % ### STATIONS TO BE DOWN-WEIGHTED ###
    elseif strcmp(str(1:min(length(str),29)),'STATIONS TO BE DOWN-WEIGHTED:')
        temp_str = textscan(str, '%s', 'CommentStyle', '#'); % parse line
        multiples = str2double(temp_str{1}{end});
        for nex = 1 : multiples
            str = fgetl(fid); % Get line
            temp_str = textscan(str, '%s', 'CommentStyle', '#'); % parse line
            switch size(temp_str{1}, 1)
                case 2 % Station name + coeff.
                    station_name_str = temp_str{1}{1};
                    coeff_str = temp_str{1}{2};
                case 3 % Station name (with blank) + coeff.
                    [flag_ok, station_name_str] = checkStatNameArrayFields(temp_str, 1, 2, str, blank_replace_char);
                    coeff_str = temp_str{1}{3};
                otherwise  % Error case
                    error('STATIONS TO BE DOWN-WEIGHTED: Invalid input argument.');
            end
            [flag_ok, station_name_str] = checkStatNameLength(station_name_str);
            ini_opt.stat_dw(nex,:) = station_name_str;
            ini_opt.stat_co(nex,:) = coeff_str;
        end


    % ### BASELINES TO BE EXCLUDED ###
    elseif strcmp(str(1:min(length(str),25)),'BASELINES TO BE EXCLUDED:')
        temp_str = textscan(str, '%s', 'CommentStyle', '#'); % parse line
        multiples = str2double(temp_str{1}{end});
        for nex = 1 : multiples
            str = fgetl(fid); % Get line
            temp_str = textscan(str, '%s', 'CommentStyle', '#'); 
            switch size(temp_str{1}, 1)
                case 2 % Two station names
                    station_name_1_str = temp_str{1}{1};
                    station_name_2_str = temp_str{1}{2};
                    
                case 3 % Two station names (one with a blank)
                    [flag_ok_1, station_name_1_str] = checkStatNameArrayFields(temp_str, 1, 2, str, blank_replace_char, false);
                    [flag_ok_2, station_name_2_str] = checkStatNameArrayFields(temp_str, 2, 3, str, blank_replace_char, false);
                    if flag_ok_1 && ~flag_ok_2
                        [flag_ok, station_name_1_str] = checkStatNameArrayFields(temp_str, 1, 2, str, blank_replace_char);
                        station_name_2_str = temp_str{1}{3};
                    elseif ~flag_ok_1 && flag_ok_2
                        [flag_ok, station_name_2_str] = checkStatNameArrayFields(temp_str, 2, 3, str, blank_replace_char);
                        station_name_1_str = temp_str{1}{1};
                    elseif flag_ok_1 && flag_ok_2
                        error('BASELINES TO BE EXCLUDED: Invalid arguments in OPT file.');
                    elseif ~flag_ok_1 && ~flag_ok_2
                        error('BASELINES TO BE EXCLUDED: Invalid arguments in OPT file.');
                    else
                        error('BASELINES TO BE EXCLUDED: Invalid arguments in OPT file.');
                    end
                    
                case 4 % Two station names (both with blanks)
                    [flag_ok, station_name_1_str] = checkStatNameArrayFields(temp_str, 1, 2, str, blank_replace_char);
                    [flag_ok, station_name_2_str] = checkStatNameArrayFields(temp_str, 3, 4, str, blank_replace_char);
                    
                otherwise  % Error case
                    error('BASELINES TO BE EXCLUDED: Invalid arguments in OPT file.');
            end
            % Check length of station names:
            [flag_ok, station_name_1_str] = checkStatNameLength(station_name_1_str);
            [flag_ok, station_name_2_str] = checkStatNameLength(station_name_2_str);
            % Save data:
            ini_opt.bas_excl(nex).sta1 = station_name_1_str;
            ini_opt.bas_excl(nex).sta2 = station_name_2_str;
            bas_excl(nex,:) = [station_name_1_str, ' ', station_name_2_str];
        end


    % ### NO CABLE CAL ###
    elseif strcmp(str(1:min(length(str),13)),'NO CABLE CAL:')
        temp_str = textscan(str, '%s', 'CommentStyle', '#'); % parse line
        multiples = str2double(temp_str{1}{end});
        for nex = 1 : multiples
            str = fgetl(fid); % Get line
            temp_str = textscan(str, '%s', 'CommentStyle', '#'); % parse line
            switch size(temp_str{1}, 1)
                case 1
                    station_name_str = temp_str{1}{1};
                case 2
                    [flag_ok, station_name_str] = checkStatNameArrayFields(temp_str, 1, 2, str, blank_replace_char);
                otherwise  % Error case
                    error('NO CABLE CAL: Invalid input argument.');
            end
            [flag_ok, station_name_str] = checkStatNameLength(station_name_str);
            ini_opt.no_cab(nex,:) = station_name_str;
        end
        
        
    % ### BASELINE DEPENDENT CLOCK OFFSET ###    
    elseif contains(str,'+BASELINE-DEPENDENT CLOCK OFFSET')
        str = fgetl(fid);
        nex=0;
        while ~contains(str,'-BASELINE-DEPENDENT CLOCK OFFSET')
            temp_str = textscan(str, '%s', 'CommentStyle', '#') ;
            if ~isempty(temp_str{1})
                nex=nex+1;
                switch size(temp_str{1}, 1)
                    case 2 % Two station names
                        station_name_1_str = temp_str{1}{1};
                        station_name_2_str = temp_str{1}{2};

                    case 3 % Two station names (one with a blank)
                        [flag_ok_1, station_name_1_str] = checkStatNameArrayFields(temp_str, 1, 2, str, blank_replace_char, false);
                        [flag_ok_2, station_name_2_str] = checkStatNameArrayFields(temp_str, 2, 3, str, blank_replace_char, false);
                        if flag_ok_1 && ~flag_ok_2
                            [flag_ok, station_name_1_str] = checkStatNameArrayFields(temp_str, 1, 2, str, blank_replace_char);
                            station_name_2_str = temp_str{1}{3};
                        elseif ~flag_ok_1 && flag_ok_2
                            [flag_ok, station_name_2_str] = checkStatNameArrayFields(temp_str, 2, 3, str, blank_replace_char);
                            station_name_1_str = temp_str{1}{1};
                        elseif flag_ok_1 && flag_ok_2
                            error('BASELINE-DEPENDENT CLOCK OFFSET: Invalid arguments in OPT file.');
                        elseif ~flag_ok_1 && ~flag_ok_2
                            error('BASELINE-DEPENDENT CLOCK OFFSET: Invalid arguments in OPT file.');
                        else
                            error('BASELINE-DEPENDENT CLOCK OFFSET: Invalid arguments in OPT file.');
                        end

                    case 4 % Two station names (both with blanks)
                        [flag_ok, station_name_1_str] = checkStatNameArrayFields(temp_str, 1, 2, str, blank_replace_char);
                        [flag_ok, station_name_2_str] = checkStatNameArrayFields(temp_str, 3, 4, str, blank_replace_char);

                    otherwise  % Error case
                        error('BASELINE-DEPENDENT CLOCK OFFSET: Invalid arguments in OPT file.');
                end

                % Check length of station names:
                [flag_ok, station_name_1_str] = checkStatNameLength(station_name_1_str);
                [flag_ok, station_name_2_str] = checkStatNameLength(station_name_2_str);
                % Save data:
                ini_opt.bdco_est(nex).sta1 = station_name_1_str;
                ini_opt.bdco_est(nex).sta2 = station_name_2_str;
             end
              str = fgetl(fid);
        end
    end
end

if remove_sprecial_stations
    ini_opt.sta_excl = [ini_opt.sta_excl;char(stations_to_be_removed(station_remove_flags))];
    
    if isfield(ini_opt, 'sta_excl_start')
        ini_opt.sta_excl_start = [ini_opt.sta_excl_start, zeros(1,sum(station_remove_flags))];
    else
        ini_opt.sta_excl_start = zeros(1,sum(station_remove_flags));
    end
end



% ##### Remove clock breaks of excluded stations #####

% Loop over all excluded stations
for i_excl_stat = 1 : size(ini_opt.sta_excl, 1)
	
    % init.:
    if exist('clk_break_tmp', 'var')
        clear clk_break_tmp
    end
    clk_break_tmp.stat_name = [];
	clk_break_tmp.mjd = [];
    i_clk_br_tmp = 0;
        
    % Loop over all clock breaks
    for i_clk_br = 1 : length(ini_opt.clk_break)
        
        % Check the station name
        if strcmp(ini_opt.clk_break(i_clk_br).stat_name, ini_opt.sta_excl(i_excl_stat,:))
            % Only, if the station is completely excluded:
            if ini_opt.sta_excl_start(i_excl_stat) == 0
                continue;
            end
        end % if strcmp(ini_opt.clk_break(i_clk_br).stat_name, ini_opt.sta_excl(i_clk_br,:))
        
        i_clk_br_tmp = i_clk_br_tmp + 1;
        clk_break_tmp(i_clk_br_tmp, :).stat_name   = ini_opt.clk_break(i_clk_br).stat_name;
        clk_break_tmp(i_clk_br_tmp, :).mjd         = ini_opt.clk_break(i_clk_br).mjd;
        
    end % for i_clk_br = 1 : length(ini_opt.clk_break)
    
    ini_opt.clk_break = clk_break_tmp;
    
end

% Set number of clock breaks:
if ~isempty([ini_opt.clk_break.mjd])
    ini_opt.num_clk_breaks = length(ini_opt.clk_break);
else
    ini_opt.num_clk_breaks = 0;
end


% ##### Close OPT file #####
fclose(fid);

end

%% ##### Sub-routines #####
function [flag_ok, stat_name_str] = checkStatNameArrayFields(temp_str, arg1_id, arg2_id, line_str, blank_replace_char, varargin)
% Checks two fields of an string-array, if they together represent one VLBI station name. 
% Only one single blank or - within a station name is permitted (leading/trailing blanks are not counted)!

    % Init:
    flag_ok = false;
    stat_name_str = '';
    
    % Optional input argument:
    switch(nargin)
       
        case 5
            flag_enable_error = true;
        case 6
            flag_enable_error = varargin{1};
        otherwise
            error('Wrong number of input arguments.');
    end
    
    % Prepare the input string "line_str" for further use:
    if min(arg1_id, arg2_id) ~= 1
        str_id = strfind(line_str, temp_str{1}{min(arg1_id, arg2_id)});
        if length(str_id) > 1
            str_id = 0;
            tmp_str_id = 1;
            k = 1;
            flag_str_id_found = false;
            while ~flag_str_id_found
                str_id = str_id + 1;
                switch line_str(str_id)

                    case ' '
                        continue;

                    case temp_str{1}{tmp_str_id}(k)
                        if tmp_str_id == min(arg1_id, arg2_id)
                            flag_str_id_found = 1;
                            % str_id found => Leave loop!
                        end
                        if k == length(temp_str{1}{tmp_str_id})
                            k = 0;
                            tmp_str_id = tmp_str_id + 1;
                        end
                        k = k + 1;

                    otherwise
                        error('Error while evaluating the string "line_str".');
                        
                end % switch
            end % while
        end % length(str_id) > 1
        
        line_str = line_str(str_id : end);
    end

    tmp_str_1 = temp_str{1}{arg1_id};
    tmp_str_2 = temp_str{1}{arg2_id};
    
    if ((length(tmp_str_1) + length(tmp_str_2)) <= 7) && (strcmp(line_str(length(tmp_str_1)+1), ' ')) && ~(strcmp(line_str(length(tmp_str_1) + 2), ' '))
        stat_name_str = [tmp_str_1, blank_replace_char, tmp_str_2];
        flag_ok = true;
    elseif ((length(tmp_str_1) + length(tmp_str_2)) <= 7) && (strcmp(line_str(length(tmp_str_1)+1), '-'))
        stat_name_str = [tmp_str_1,'-', tmp_str_2];
        flag_ok = true;
    else
        if flag_enable_error
            error('readOPT: Invalid station name in OPT file.');
        end
    end
	
end

function [flag_ok, stat_name_str] = checkStatNameLength(stat_name_str)
% Checks if the station name is an 8 char. sring. 
% If < 8 char. => Add blanks!
% If > 8 char. => Error!

    % Init.:
    flag_ok = true;

    % Check max length (8 char.):
    if length(stat_name_str) > 8
        flag_ok = false;
        error('readOPT: Invalid input argument (station name with >8 char.).'); 
    end
    
    % Add blanks if required to get 8 char.:
    stat_name_str = sprintf('%-*s' , 8 , stat_name_str); % 8 char station name => add blanks, if necessary

end


% TO Do:
% - 1.) Check "excl. baselines" => case 3; similar station names (hoe to
% choose sub-str for "line_str"?
% - 2.) Write warning to CW, WHENEVER a blank is found in an OPT file??? =>
% Correct this in the according OPT file manually! ??
