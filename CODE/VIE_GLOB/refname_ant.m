% ************************************************************************
%   Description:
%   Loop over all sessions to get all names of antennas and their intervals
%
%
%   Input:
%      path              path to LEVEL2 data
%      ses               names of sessions
%      antbr             antenna breaks from an external file (provided by IVS analysis coordinator (A.Nothnagel, Bonn))
%                        'VLBI-DISCONT.txt' converted to a mat.file
%     maxRMS             maximal RMS for this adjustment
%     dir_in             name of the input directory
%
%   Output:                
%     refname            names of all antennas
%     refantbr           array with all breaks for all stations
%     antactiv           matrix with 1/0; columns are sessions; rows are
%                        antennas. 1 - antenna observed in this session, 0 - antenna didn't take
%                        part in this session
%     mjd_all            vector with mjd of all sessions
%       
%   Coded for VieVS: 
%   21 Jul 2011 by Hana Spicakova
%
%   Revision: 
%   16 Sep 2011 by Hana Spicakova    maxRMS and dir_in added as input
%   parameters
%%



function [refname,refantbr,antactiv,mjd_all]=refname_ant(path, ses, antbr, maxRMS, dir_in)

interv=[];
mjd_all=[];
refname=[''];

lse=size(ses,2);

for ise=1:lse % loop over all sessions
    load ([path ses{ise} '_an_glob.mat']);

    antenna = glob1.an;  
           
    % if there is a new station, its name and coordinates are added
    asize = length(antenna.x);  % number of antennas in this session
    anames=antenna.name;
    mjd=antenna.firstscan_mjd;
    for i=1:asize  % antennas in one session
        aname=anames(i,:);
        fa=strcmp(cellstr(aname),cellstr(refname)); % find the antenna
        nra=find(fa,1); % number of the antenna
        fa=sum(fa); % 1-found / 0-not found
        ln = size(refname,1); % final number of stations (breaks are not considered)
      
        for j=1:length(antbr)
            r(j)= strcmp(cellstr(aname),cellstr(antbr(j).name));
        end
        j=find(r==1);
        if fa==0
            refantbr(ln+1,:).name = aname;
            % breaks for the stations from an external file 'VLBI-DISCONT.txt'
            if isempty(j)==0
                refantbr(ln+1,:).break = antbr(j).break;
                refantbr(ln+1,:).break(length(refantbr(ln+1,:).break)+1) = 99999;
                refantbr(ln+1,:).interv = zeros(1,length(refantbr(ln+1,:).break)-1);
                if strcmp(refantbr(ln+1,:).name,'GILCREEK')
                    refantbr(ln+1,:).break(length(refantbr(ln+1,:).break)) = 54100; % GILCREEK stopped observations on this day
                end
            else
                refantbr(ln+1,:).break = [0 99999];
                refantbr(ln+1,:).interv = 0;
            end
            
            % active station for a plot
            antactiv(ln+1,ise)=1;
        else
            antactiv(nra,ise)=1;
        end
                
        if fa~=0; nrant=nra; else nrant=ln+1; end
        for k=1:length(refantbr(nrant).break)-1
            gg = find(refantbr(nrant).break(k)<mjd && mjd<refantbr(nrant).break(k+1));
            if isempty(gg); gg=0; end
            interv(k)=gg;
        end

        % number of sessions with this station in the intervales defined in
        % 'VLBI-DISCONT.txt'
        refantbr(nrant,:).interv = refantbr(nrant,:).interv+interv;
        clear interv
       
        ln = length(refantbr); % final number of stations (breaks are not considered)
        for j = 1: ln
            refname(j,:) = refantbr(j,:).name;
        end
    end
    mjd_all=[mjd_all mjd];
    clear glob1
end

isln = exist('ln','var');
if isln==0
    fprintf('\n\n No sessions with RMS < %5.2f found in the directory %s! \n',maxRMS, dir_in);
end


