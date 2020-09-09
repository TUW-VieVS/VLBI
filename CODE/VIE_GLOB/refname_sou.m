% ************************************************************************
%   Description:
%   Loop over all sessions to get all names of sources
%
%
%   Input:
%      path              path to LEVEL2 data
%      ses               names of sessions
%      fixedsou          names of sources, which will be fixed to apriori
%                        coordinates (= deleted from N matrix)
%      reducsou          names od sources, which will be session-wise
%                        reduced
%
%
%   Output:     
%     qrefname           names of sources going to global adjustment
%     souactiv           matrix with numobs/0; columns are sessions; rows are
%                        sources going to global adjustment. If the source was 
%                        observed in this session numobs is saved, 0 - source didn't take
%                        part in this session
%       
%   Coded for VieVS: 
%   21 Jul 2011 by Hana Spicakova
%
%   Revision: 
%   10 Oct 2012 by Hana Krasna   if a source is observed in a session,
%       instead of 1 the number of observations is stored in the souactiv
%       matrix (needed for the official CRF catalogue output)
%   24 Mar 2013 by Hana Krasna   IERS name is stored (instead of the name in NGS which could be also IVS/common name)
%	10 Apr 2017 by David Mayer IVS names are also saved
%%




function [qrefname, souactiv] = refname_sou(path, ses, fixedsou, reducsou)

qrefname.IERS=[''];
qrefname.IVS=[''];


s=0;
lse=size(ses,2);

for ise=1:lse
    load ([path ses{ise} '_par_glob.mat'])

    qsize=length(glob2.opt.source);
    for iso=1:qsize
        qnameIERS = glob2.opt.source(iso).IERSname;
        qnameIVS  = glob2.opt.source(iso).name;
        nrq=find(strcmp(cellstr(qnameIERS),cellstr(qrefname.IERS)),1); % find the source
        if isempty(nrq)
            s=s+1;
            qrefname.IERS(s,:)=qnameIERS;
            qrefname.IVS(s,:)=qnameIVS;
            souactiv(s,ise)=glob2.opt.source(iso).total_obs;
        else
            souactiv(nrq,ise)=glob2.opt.source(iso).total_obs;
        end
    end

clear nrq
end

%qrefname_all=qrefname;
%souactiv_all=souactiv;


% delete the fixed sources from qrefname and souactiv
nrq=[];
l=1;
for i=1:size(fixedsou,1)
    fqname=fixedsou(i,:);
    j=find(strcmp(cellstr(fqname),cellstr(qrefname.IERS)),1);
    if ~isempty(j)
       nrq(l)=j;
       l=l+1;
    end
    clear k j nrq1
end

qrefname.IERS(nrq,:)=[];
qrefname.IVS(nrq,:)=[];
souactiv(nrq,:)=[];
    
clear nrq
   
% delete the session-wise reduced sources from qrefname and souactiv
nrq=[];
l=1;
for i=1:size(reducsou,1)
    fqname=reducsou(i,:);
    j=find(strcmp(cellstr(fqname),cellstr(qrefname.IERS)),1);
    if ~isempty(j)
        nrq(l)=j;
        l=l+1;
    end
    clear k j nrq1
end

qrefname.IERS(nrq,:)=[];
qrefname.IVS(nrq,:)=[];

souactiv(nrq,:)=[];






