% ************************************************************************
% Description:
%    Tool to make process-lists for VieVS. Reads information about the
%    sessions from masterfiles in $VIEVSROOT/DATA/MASTER.
%    Just process sessions which are in the masterfiles and in the folder
%    /DATA/NGS/YEAR/
%
%  Examples:
%   To get all R1 & R4 sessions run:
%      [list,sessionnames]=mk_list('R1','R4')
%   To get all R1 & R4 sessions run from the period 2005-2010, excluding
%   stations given in '../DATA/MASTER/exclude.txt':
%      [list,sessionnames]=mk_list('R1','R4','YEARS',2005:2010,...
%             'EXCLUDE','../DATA/MASTER/exclude.txt')
%   To get all EURO experiments where ONSALA60 has participated run:
%      [list,sessionnames]=mk_list('EURO','REQSTAT','On')
%
% Usage:
%    [list,sessionnames]=mk_list(s1,s1,...)
%      where s1,s2,... are the names of the types of sessions you want to
%      include (R1, R4, EURO, etc.). Use all for all sessions
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
%    non-intensiv sessions
%    20 Apr 2015 by David Mayer: MINSTANUM added
%    31 MAR 2017 by David Mayer: corrected NOINT option 
%
% ************************************************************************

function [list,names]=mk_list(varargin)


years=1979:(clock*[1;0;0;0;0;0]);
excl={};
expes={};
reqstat={};
allv=0;
nointv=0;
minstanum=0;
a=1;
while a<=nargin
    str=varargin{a};
    switch upper(str)
        case 'YEARS'
            a=a+1;
            years=varargin{a};
        case {'EXCLUDE','EXCL'}
            a=a+1;
            excl=[excl varargin{a}];
        case {'ALLVERS','ALLVERSIONS'}
            allv=1;
        case {'REQSTAT','ONLY'}
            a=a+1;
            reqstat=[reqstat varargin{a}];
        case {'NOINT'}
            nointv=1;
        case 'MINSTANUM'
            a=a+1;
            minstanum=varargin{a};
        otherwise
            expes=[expes str];
            switch str
                case 'EUR'
                    expes=[expes 'EURO'];
                case 'EURO'
                    expes=[expes 'EUR'];
            end
    end
    a=a+1;
end

excllist=[];
for a=1:length(excl)
    fid=fopen(excl{a},'r');
    while ~feof(fid)
        str=[fgetl(fid) '           '];
        if (str(1)=='$')&(str(2)~='$')
            tmp=str(2:10);
            tmp(find(tmp==' '))='_';
            excllist=[excllist;tmp];
        end
    end
    fclose(fid);
end

names=[];
list=[];
for a=1:length(years)
    year=years(a);
    yrsh=year-100*floor(year/100);
    pt=sprintf('../DATA/NGS/%4d/',year);
    for intf=0:1
        if intf
            if nointv==1
                continue
            end
            mfil=sprintf('../DATA/MASTER/master%02d-int.txt',yrsh);
        else
            mfil=sprintf('../DATA/MASTER/master%02d.txt',yrsh);
        end
        if exist(mfil)
            fid=fopen(mfil,'r');
            while ~feof(fid)
                str=[fgetl(fid) '   '];
                
                
                if str(1)=='|'
                    id=find(str=='|');
                    
                    expnam=str(id(2)+1:id(3)-1);
                    incl=0;
                    for a=1:length(expes)
                        ids=id(2)+1:id(2)+length(expes{a});
                        if strcmp(str(ids),expes{a})&...
                                (length(str2num(str(ids(end)+1)))==1)
                            incl=1;
                            break
                        end
                        if strcmp(upper(expes{a}),'ALL')
                            incl=1;
                        end
                    end
                    date=str(id(3)+1:id(4)-1);
                    date(find(date==' '))=[];
                    code=str(id(12)+1:id(13)-1);
                    code(find(code==' '))=[];
                    if length(code)==1
                        code=[code '_'];
                    elseif length(code)==0
                        code='__';
                    end
                    sess=sprintf('%02d%s%s',yrsh,date,code);
                    a=1;
                    while  (a<=size(excllist,1))&incl
                        if strcmp(excllist(a,:),sess)
                            incl=0;
                        end
                        a=a+1;
                    end
                    if incl&(length(reqstat)>0)
                        tmp=str(id(7)+1:id(8));
                        ids=find(tmp==' ');
                        tmp=tmp(1:ids(1)-1);
                        for a=1:length(reqstat)
                            if isempty(findstr(reqstat{a},tmp))
                                incl=0;
                            end
                        end
                        
                    end
                    if incl&(minstanum>0)
                        tmp=deblank(str(id(7)+1:id(8)-1)); % remove blanks from string
                        find_dash = strfind(tmp, ' -');
                        if ~isempty(find_dash)
                            tmp = tmp(1:find_dash-1); % delete 'removed stations' from the string
                        end
                        if length(tmp)/2 < minstanum
                            incl=0;
                        end
                    end
                    if incl
                        fils=dir([pt sess '*']);
                        listtmp=[];
                        ver=[];
                        for b=1:length(fils)
                            listtmp=[listtmp;sprintf('%4d/%s',year,fils(b).name)];
                            ver=[ver str2num(fils(b).name(12:14))];
                        end
                        if allv
                            list=[list;listtmp]
                            for b=1:size(listtmp,1);
                                names=[names;expnam];
                            end
                        elseif ~isempty(ver)
                            [m,mid]=max(ver);
                            list=[list;listtmp(mid,:)];
                            names=[names;expnam];
                        end
                    end
                end
            end
            fclose(fid);
        end
    end
end
