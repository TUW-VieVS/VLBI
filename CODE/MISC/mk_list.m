% ************************************************************************
% Description:
%    Tool to make process-lists for VieVS. Reads information about the
%    sessions from masterfiles in $VIEVSROOT/DATA/MASTER.
%    Just process sessions which are in the masterfiles and in the folder
%    /DATA/NGS/YEAR/.
%    Attention: mk_list.m works for NGS-files only, vgosDB is not
%    supported!
%
%  Examples:
%   To get all R1 & R4 sessions run:
%      [list,sessionnames]=mk_list('R1','R4')
%   To get all R1 & R4 sessions run from the period 2005-2010, excluding
%   stations given in '../DATA/MASTER/exclude.txt':
%      [list,sessionnames]=mk_list('R1','R4','YEARS',2005:2010,...
%             'EXCLUDE','../DATA/MASTER/exclude.txt')
%   To get all EURO experiments where ONSALA60 has participated in:
%      [list,sessionnames]=mk_list('EURO','REQSTAT','On')
%
% Usage:
%    [list,sessionnames]=mk_list(s1,s1,...)
%      where s1,s2,... are the names of the types of sessions you want to
%      include (R1, R4, EURO, etc.). Use all for all sessions
%    [list,sessionnames]=mk_list(...'VGOSDB') or process_list=mk_list(...,'NGS')
%      to define whether vgosDB or NGS input files are used
%    [list,sessionnames]=mk_list(...,'YEARS',yrs)
%      Only use sessions from the years 'yrs'.
%    [list,sessionnames]=mk_list(...,'EXCLUDE',excllist)
%      Exclude sessions listed in the file 'excllist'. The format of this
%      list should be similar to the list found at:
%            http://vlbi.geod.uni-bonn.de/IVS-AC/data/exclude.txt
%    [list,sessionnames]=mk_list(...,'ALLVERSIONS')
%       Include all availible versions of the NGS-files in the list, for
%       the cases several versions are availible (defult is only to include
%       the latest version).
%    [list,sessionnames]=mk_list(...,'REQSTAT',{stats})
%       Only include sessions where all of the stations 'stats' are
%       participating. The station names givben in 'stats' should be the
%       2-letter ns-codes defined in the file:
%          ftp://cddis.gsfc.nasa.gov/pub/vlbi/ivscontrol/ns-codes.txt
%    [list,sessionnames]=mk_list(...,'NOINT')
%       Include all availible sessions of the NGS-files, which are
%       non-intensive sessions
%    [list,sessionnames]=mk_list(...,'MINSTANUM',4)
%       Include all available sessions which have 4 or more stations
%
% Output:
%     list      (n*19 string)  List of NGS-files
%     sessnames (n*6 string)   6-letter naames of the sessions
%
% Coded for VieVS:
%    12 Feb 2010 by Tobias Nilsson
%
% Revision:
%    19 Jan 2015 by Caroline Schönberger: KEYWORD 'NOINT' inlcude all
%    non-intensive sessions
%    20 Apr 2015 by David Mayer: MINSTANUM added
%    31 Mar 2017 by David Mayer: corrected NOINT option 
%    23 Aug 2018 by Daniel Landskron: option ONLYINT added
%    22 Oct 2018 by Daniel Landskron: also enabled for vgosDB
%    07 Dec 2018 by Daniel Landskron: command window output extended and bug corrected
%
% ************************************************************************

function [list,names]=mk_list(varargin)


% set initial values 
session_data_format = 'VGOSDB';
years = 1979:(clock*[1;0;0;0;0;0]);
exclude = {};
expes ={};
reqstat = {};
all_versions = 0;
no_intensives = 0;
only_intensives = 0;
min_station_num = 0;


i_input = 1;
while i_input <= nargin
    
    str=varargin{i_input};
    switch upper(str)
        case 'NGS'
            session_data_format = 'NGS';
        case 'YEARS'
            i_input = i_input+1;
            years = varargin{i_input};
        case {'EXCLUDE','EXCL'}
            i_input = i_input+1;
            exclude = [exclude varargin{i_input}];
        case {'ALLVERS','ALLVERSIONS'}
            all_versions = 1;
        case {'REQSTAT','ONLY'}
            i_input = i_input+1;
            reqstat = [reqstat varargin{i_input}];
        case {'NOINT'}
            no_intensives = 1;
        case {'ONLYINT'}
            only_intensives = 1;
        case 'MINSTANUM'
            i_input = i_input+1;
            min_station_num = varargin{i_input};
        otherwise
            expes = [expes str];
            switch str
                case 'EUR'
                    expes = [expes 'EURO'];
                case 'EURO'
                    expes = [expes 'EUR'];
            end
    end
    i_input = i_input+1;
    
end

excllist = [];
for a = 1:length(exclude)
    fid = fopen(exclude{a},'r');
    while ~feof(fid)
        str = [fgetl(fid) '           '];
        if (str(1)=='$')  &&  (str(2)~='$')
            tmp = str(2:10);
            tmp(tmp==' ') = '_';
            excllist = [excllist;tmp];
        end
    end
    fclose(fid);
end

names = [];
list = [];
for i_year = 1:length(years)
    year = years(i_year);
    year_short = year-100*floor(year/100);
    
    if strcmpi(session_data_format,'VGOSDB')
        path_session = sprintf('../DATA/vgosDB/%4d/',year);
    else
        path_session = sprintf('../DATA/NGS/%4d/',year);
    end
        
    for intf = 0:1
        if intf
            if no_intensives==1
                continue
            end
            path_master_file = sprintf('../DATA/MASTER/master%02d-int.txt',year_short);
        else
            if only_intensives==1
                continue
            end
            path_master_file = sprintf('../DATA/MASTER/master%02d.txt',year_short);
        end
        if exist(path_master_file,'file')
            
            fid = fopen(path_master_file,'r');
            while ~feof(fid)
                str = [fgetl(fid) '   '];
                
                
                if str(1)=='|'
                    id = find(str=='|');
                    
                    expnam = str(id(2)+1:id(3)-1);
                    incl = 0;
                    for i_expes = 1:length(expes)
                        ids = id(2)+1:id(2)+length(expes{i_expes});
                        if strcmp(str(ids),expes{i_expes})  &&  (length(str2num(str(ids(end)+1)))==1)
                            incl = 1;
                            break
                        end
                        if strcmpi(expes{i_expes},'ALL')
                            incl = 1;
                        end
                    end
                    date = str(id(3)+1:id(4)-1);
                    date(date==' ') = [];
                    code = str(id(12)+1:id(13)-1);
                    code(code==' ') = [];
                    if length(code)==1
                        code = [code '_'];
                    elseif isempty(code)
                        code = '__';
                    end
                    sess = sprintf('%02d%s%s',year_short,date,code);
                    
                    a = 1;
                    while  (a<=size(excllist,1))  &&  incl
                        if strcmp(excllist(a,:),sess)
                            incl = 0;
                        end
                        a = a+1;
                    end
                    
                    if incl  &&  (~isempty(reqstat))
                        tmp = str(id(7)+1:id(8));
                        ids = find(tmp==' ');
                        tmp = tmp(1:ids(1)-1);
                        for i_reqstat = 1:length(reqstat)
                            if isempty(findstr(reqstat{i_reqstat},tmp))
                                incl = 0;
                            end
                        end
                        
                    end
                    if incl  &&  (min_station_num>0)
                        tmp = deblank(str(id(7)+1:id(8)-1)); % remove blanks from string
                        find_dash = strfind(tmp, ' -');
                        if ~isempty(find_dash)
                            tmp = tmp(1:find_dash-1); % delete 'removed stations' from the string
                        end
                        if length(tmp)/2 < min_station_num
                            incl = 0;
                        end
                    end
                    if incl
                        if strcmpi(session_data_format,'VGOSDB')
                           % unpack vgosDB *.tgz or *.tar.gz
                           sess_tgz = [sess,'.tgz'];
                           sess_targz = [sess,'.tar.gz'];
                           if exist([path_session sess_tgz],'file')
                              untar([path_session sess_tgz],path_session);
                           elseif exist([path_session sess_targz],'file')
                              untar([path_session sess_targz],path_session);
                           else
                           fprintf('ERROR: %s is not available in /DATA/vgosDB!\n',sess_tgz);
                           end 
                         end
                        
                        session_file = dir([path_session sess '*']);   % can be >1, in case there are several versions
                        listtmp = [];
                        ver=[];
                        for i_session_file = 1:length(session_file)
                            if length(session_file(i_session_file).name) <= 9   % because there are some sessions with '_p1' or '.tmp', which shall not be considered
                                listtmp = [listtmp;sprintf('%4d/%s',year,session_file(i_session_file).name)];
                            end
                            if strcmpi(session_data_format,'NGS')
                                listtmp = [listtmp;sprintf('%4d/%s',year,session_file(i_session_file).name)];
                                ver = [ver str2num(session_file(i_session_file).name(12:14))];
                            elseif strcmpi(session_data_format,'VGOSDB')    
                                all_versions = 0;   % for vgosDB file only the latest version is used
                            end
                        end
                        if all_versions
                            list = [list;listtmp];
                            for b = 1:size(listtmp,1)
                                names = [names;expnam];
                            end
                        elseif ~isempty(ver)
                            [~,i_max] = max(ver);
                            list = [list;listtmp(i_max,:)];
                            names = [names;expnam];
                        elseif strcmpi(session_data_format,'VGOSDB')    && ~isempty(listtmp)
                            list = [list;listtmp ' [vgosDB]'];
                            names = [names;expnam];
                        end
                        % remove the unpacked vgosDB folder
                        if strcmpi(session_data_format,'VGOSDB') && exist([path_session sess],'dir')
                        rmdir([path_session sess], 's');
                        end
                    end
                end
                
            end
            fclose(fid);
            
        else
            
            % give a message that the necessary master file can't be found, but make exceptions:
            % - there are no Intensive master files before 1992
            if ~contains(path_master_file,'-int')
                fprintf('%s%s%s\n','The Intensives master file ',path_master_file,' is missing! However, If you don''t consider Intensive sessions, this is no problem.');
            elseif contains(path_master_file,'-int') && years(i_year)>=1992 
                fprintf('%s%s%s\n','The master file ',path_master_file,' is missing! However, If you only consider Intensive sessions, this is no problem.');
            end
            
        end
    end
end

fclose('all');
