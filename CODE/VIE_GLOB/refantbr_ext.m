% ************************************************************************
%   Description:
%   Find the apriori coordinates for stations for all intervals
%
%
%   Input:
%      refantbr          array with all breaks for stations in global adjustment
%      refname           names of all antennas
%      trf               trf catalogue
%      trffil            specification of the trf catalogue
%      path              path to LEVEL2 data
%      ses               names of sessions
%
%   Output:                
%      X0                vector with x apriori coordinates of stations in [m]
%      Y0                vector with y apriori coordinates of stations in [m]
%      Z0                vector with z apriori coordinates of stations in [m]
%      brstart           vector with time when the coordinate interval starts [mjd]
%      brend             vector with time when the coordinate interval ends [mjd]
%      refantbr          array with all breaks for stations in global
%                        adjustment and the apriori coordinates
%
%   Coded for VieVS: 
%   21 Jul 2011 by Hana Spicakova
%
%   Revision: 
%   07 Aug 2012 by Hana Krásná: change to be compatible with the
%                       superstation file
%   04 Oct 2012 by Hana Krásná: small bug fixed which occured if there
%                  were no observations in the first interval of station
%                  coordinates
%   22 Nov 2012 by Hana Krásná: small bug fixed which caused the program to
%            crashed if using the backup vievs-coordinates
%   25 Mar 2013 by Hana Krásná: epochs added to refantbr
%   11 Dec 2015 by Andreas Hellerschmied: Bug fix (cellstr)
%%

function [X0,Y0,Z0,brstart,brend,refantbr] = refantbr_ext(refantbr,refname,trf,trffil,path,ses)

 
 
    for i=1:size(refname,1)
        clear j r aname 
        aname=refname(i,:);
        
        idblank=find(aname(1:max(find(aname(:)~=' ')))==' ');
        aname_=aname;
        aname_(idblank)='_';
        
        for k = 1:length(trf)
            r(k)= strcmp(aname_,trf(k).name);
        end
        j=find(r==1);
       
        
        % save the antenna apriori coordinates for all breaks from trf
        % catalogue
        if isempty(j)==0
            if ~isempty(trf(j).(trffil{2})) % if it is in chosen catalogue
                for k = 1:length(trf(j).(trffil{2}).break)
                    if ~isempty(trf(j).(trffil{2}).break(k).start) % HLOUPY - kvuli Dss15 . vtrf2008 se musi predelat! HANA
                        refantbr(i,:).break_apr(:,k) = trf(j).(trffil{2}).break(k).start;
                        refantbr(i,:).x_apr(:,k) = trf(j).(trffil{2}).break(k).x;
                        refantbr(i,:).y_apr(:,k) = trf(j).(trffil{2}).break(k).y;
                        refantbr(i,:).z_apr(:,k) = trf(j).(trffil{2}).break(k).z;
                        refantbr(i,:).vx_apr(:,k) = trf(j).(trffil{2}).break(k).vx;
                        refantbr(i,:).vy_apr(:,k) = trf(j).(trffil{2}).break(k).vy;
                        refantbr(i,:).vz_apr(:,k) = trf(j).(trffil{2}).break(k).vz;
                        refantbr(i,:).epoch(:,k) = trf(j).(trffil{2}).break(k).epoch;
                    else
                        refantbr(i,:).break_apr(:,k) = trf(j).vievsTrf.break(k).start;
                        refantbr(i,:).x_apr(:,k) = trf(j).vievsTrf.break(k).x;
                        refantbr(i,:).y_apr(:,k) = trf(j).vievsTrf.break(k).y;
                        refantbr(i,:).z_apr(:,k) = trf(j).vievsTrf.break(k).z;
                        refantbr(i,:).vx_apr(:,k) = trf(j).vievsTrf.break(k).vx;
                        refantbr(i,:).vy_apr(:,k) = trf(j).vievsTrf.break(k).vy;
                        refantbr(i,:).vz_apr(:,k) = trf(j).vievsTrf.break(k).vz;
                        refantbr(i,:).epoch(:,k) = trf(j).vievsTrf.break(k).epoch;
                    end
                end
                refantbr(i,:).break_apr(:,k+1) = trf(j).(trffil{2}).break(k).end;
            elseif ~isempty(trf(j).vievsTrf) % backup vievsTRF catalogue
                for k = 1:length(trf(j).vievsTrf.break)
                    refantbr(i,:).x_apr(:,k) = trf(j).vievsTrf.break(k).x;
                    refantbr(i,:).y_apr(:,k) = trf(j).vievsTrf.break(k).y;
                    refantbr(i,:).z_apr(:,k) = trf(j).vievsTrf.break(k).z;
                    
                    if isfield(trf(j).vievsTrf.break(k), 'start')
                        refantbr(i,:).break_apr(:,k) = trf(j).vievsTrf.break(k).start;
                        refantbr(i,:).vx_apr(:,k) = trf(j).vievsTrf.break(k).vx;
                        refantbr(i,:).vy_apr(:,k) = trf(j).vievsTrf.break(k).vy;
                        refantbr(i,:).vz_apr(:,k) = trf(j).vievsTrf.break(k).vz;
                        refantbr(i,:).epoch(:,k) = trf(j).vievsTrf.break(k).epoch;
                    else
                        refantbr(i,:).break_apr(:,k) = 0;
                        refantbr(i,:).vx_apr(:,k) = 0;
                        refantbr(i,:).vy_apr(:,k) = 0;
                        refantbr(i,:).vz_apr(:,k) = 0;
                        refantbr(i,:).epoch(:,k) = 0;
                    end
                    
                    
                end
                refantbr(i,:).break_apr(:,k+1) = 99999;
            else % if it is nowhere
                refantbr(i,:).break_apr = [0 99999];
                refantbr(i,:).x_apr = 0;
                refantbr(i,:).y_apr = 0;
                refantbr(i,:).z_apr = 0;
                refantbr(i,:).vx_apr = 0;
                refantbr(i,:).vy_apr = 0;
                refantbr(i,:).vz_apr = 0;
                refantbr(i,:).epoch = 0;
            end
            
        else % if the station is not in the catalogue
            refantbr(i,:).break_apr = [0 99999];
            refantbr(i,:).x_apr = 0;
            refantbr(i,:).y_apr = 0;
            refantbr(i,:).z_apr = 0;
            refantbr(i,:).vx_apr = 0;
            refantbr(i,:).vy_apr = 0;
            refantbr(i,:).vz_apr = 0;
            refantbr(i,:).epoch = 0;
        end
        refantbr(i,:).interv_apr = zeros(1,length(refantbr(i,:).break_apr)-1);
    end
    
    % apriori intervals --> find out, how many observations were in the
    % intervals of the trf catalogue
    lse=size(ses,2);
    for ise=1:lse
        
        load ([path ses{ise} '_an_glob.mat']);

        antenna = glob1.an;
        asize = length(antenna.x);  % number of antennas in this session
        anames=antenna.name;
        mjd=antenna.firstscan_mjd;
        for i=1:asize  % antennas in one session
            aname=anames(i,:);
            fa=strcmp(cellstr(aname),cellstr(refname)); % find the antenna
            if sum(fa)~=0
                nra=find(fa,1); % number of the antenna
                for k=1:length(refantbr(nra).break_apr)-1
                    gg = find(refantbr(nra).break_apr(k)<mjd && mjd<refantbr(nra).break_apr(k+1));
                    if isempty(gg); gg=0; end
                    interv_apr(k)=gg;
                end
                refantbr(nra,:).interv_apr = refantbr(nra,:).interv_apr + interv_apr;

                % if the antenna was not in the trf catalogue, take the
                % coordinates from the antenna structure 
                if refantbr(nra,:).x_apr==0 
                    refantbr(nra,:).x_apr=antenna.x(i);
                    refantbr(nra,:).y_apr=antenna.y(i);
                    refantbr(nra,:).z_apr=antenna.z(i);
                    refantbr(nra,:).epoch=antenna.epoch(i);
                end  
            end
            clear interv_apr
        end
    end
    
    % check if number of observation in refantbr.interv is the same as in
    % refantbr.interv_apr; It is needed, if new coordinates after
    % earthquake are still not included in trf.
    for i=1:length(refantbr)
        if sum(refantbr(i).interv) ~= sum(refantbr(i).interv_apr)
            for ise=1:lse
                load ([path ses{ise} '_an_glob.mat']);
                antenna = glob1.an;
                if antenna.firstscan_mjd>refantbr(i).break_apr(end)
                    f=find(strcmp(cellstr(refantbr(i).name),cellstr(glob1.an.name))==1);
                    if ~isempty(f)
                        refantbr(i).break_apr(end+1)=99999;
                        refantbr(i).x_apr(end+1)=antenna.x(f);
                        refantbr(i).y_apr(end+1)=antenna.y(f);
                        refantbr(i).z_apr(end+1)=antenna.z(f);
                        refantbr(i).vx_apr(end+1)=antenna.vx(f);
                        refantbr(i).vy_apr(end+1)=antenna.vy(f);
                        refantbr(i).vz_apr(end+1)=antenna.vz(f);
                        refantbr(i).epoch(end+1)=antenna.epoch(f);
                        refantbr(i).interv_apr(end+1)=sum(refantbr(i).interv)-sum(refantbr(i).interv_apr);
                        fprintf('\n\n For station %s were added apriori coordinates since %5.0f, which were not in the trf catalogue! \n You can check the variable refantbr(%1.0f). \n', refantbr(i).name,refantbr(i).break_apr(end-1),i);
                        break
                    end
                end
            end
            
            for ise=1:lse
                load ([path ses{ise} '_an_glob.mat']);
                antenna = glob1.an;

                if ~isempty(find(strcmp(cellstr(refantbr(i).name),cellstr(glob1.an.name))==1))
                % if the first intervals are missing
                    if antenna.firstscan_mjd<refantbr(i).break_apr(1)

                        f=find(strcmp(cellstr(refantbr(i).name),cellstr(glob1.an.name))==1);
                        refantbr(i).break_apr=[0 refantbr(i).break_apr] ;
                        refantbr(i).x_apr = [antenna.x(f)  refantbr(i).x_apr];
                        refantbr(i).y_apr = [antenna.y(f) refantbr(i).y_apr];
                        refantbr(i).z_apr = [antenna.z(f) refantbr(i).z_apr];
                        refantbr(i).vx_apr = [antenna.vx(f) refantbr(i).vx_apr];
                        refantbr(i).vy_apr = [antenna.vy(f)  refantbr(i).vy_apr];
                        refantbr(i).vz_apr = [antenna.vz(f) refantbr(i).vz_apr];
                        refantbr(i).epoch = [antenna.epoch(f) refantbr(i).epoch];
                        refantbr(i).interv_apr = [sum(refantbr(i).interv)-sum(refantbr(i).interv_apr) refantbr(i).interv_apr];

                        fprintf('\n\n For station %s were added apriori coordinates before %5.0f, which were not in the trf catalogue! \n You can check the variable refantbr(%1.0f). \n', refantbr(i).name,refantbr(i).break_apr(2),i);
                        break
                    end
                end
            end
        end
    end
    % apriori coordinates for the new intervals (for datum definition and new TRF)
    k=0;
    for i=1:size(refname,1)
        interv_N=find(refantbr(i,:).interv>0); % index of the interval (breaks from external file)
        startbr_N=refantbr(i,:).break(interv_N); % when start these intervals
        endbr_N=refantbr(i,:).break(interv_N+1); % when end    -||-
        for j=1:length(interv_N)
            k=k+1;
            XYZinda=find(startbr_N(j)>=refantbr(i,:).break_apr); % because of a possible new break in the external file, it must be used >=
            if isempty(XYZinda) % if the date of breaks from external file is before beginning of the catalogue
                XYZind=1;
            else
                XYZind=XYZinda(end);
            end
            X0(k)=refantbr(i,:).x_apr(XYZind);
            Y0(k)=refantbr(i,:).y_apr(XYZind);
            Z0(k)=refantbr(i,:).z_apr(XYZind);
            brstart(k)=startbr_N(j);
            brend(k)=endbr_N(j);
        end
    end
