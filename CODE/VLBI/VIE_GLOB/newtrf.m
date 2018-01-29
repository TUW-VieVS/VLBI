% ************************************************************************
%   Description:
%   This function creates a .txt file with new station coordinates and 
%   velocities. This file can be directly used by VieVS, as your own TRF catalogue
%   (copy this file into VieVS/TRF/)
%
%   Input:										
%       refantbr            structure with information about
%                           discontinuities in station positions
%       globsol             structure with estimates
%       paths               paths to directories
%                           
%   Output:                
%      'trf_*.txt'
%
%   External calls: 	
%      -               					    											
%       
%   Coded for VieVS: 
%   03 Aug 2010 by Hana Spicakova
%
%   Revision: 
%   12 Jan 2011 by Hana Spicakova:  a bug corrected which occured if a
%      break in external file is before beginning of the used catalogue
%   03 Mar 2011 by Hana Spicakova: newly written
%   14 Apr 2011 by Hana Spicakova: some minor changes
%   17 Oct 2012 by Hana Krásná: Output coded in a format of an offical TRF
%               catalogue
%   21 Jun 2013 by Hana Krasna: VieVS catalogue is sorted alphabetically
%   11 Dec 2015 by Andreas Hellerschmied: Bug fix (strfind)
%   07 Jun 2017 by Hana Krasna: bug fixed by computing the new position (if the estimated catalogue has different epoch as the a priori)
%%

 function newtrf(refantbr,globsol,paths,trf)

if exist([paths.path_out 'TRF/' paths.out])~=7
	mkdir([paths.path_out 'TRF/' paths.out])
end


nst=size(globsol.antenna.refname_coord,1);

if globsol.antenna.id_dxyz == 0
   globsol.antenna.dxyz = zeros(nst,3);
end
if globsol.antenna.id_dvxyz == 0
   globsol.antenna.dvxyz = zeros(nst,3);
end


% get current date and time
curDate=clock;          % current date and time



%% Create an ASCII TRF Catalogue in an "official" format
% only if coordinates and velocities were estimated
if globsol.antenna.id_dvxyz == 1
    fidOffic=fopen([paths.path_out 'TRF/' paths.out '/trf_catalogue_' paths.L2 '.txt'],'wt');


    timeepoch=mjd2yydoysecod(globsol.antenna.epoch_break(1,1));
    timeep = timeepoch(:,1)+timeepoch(2)/365;

    fprintf(fidOffic, 'Created on %02.0f.%02.0f.%04.0f at %02.0f:%02.0f:%02.0f local time\n', curDate(3), curDate(2), curDate(1), curDate(4), curDate(5), curDate(6));
    fprintf(fidOffic,'                                  STATION POSITIONS AT EPOCH %5.1f AND VELOCITIES \n',timeep);
    fprintf(fidOffic,'                                                VLBI STATIONS \n\n\n');
    fprintf(fidOffic,'DOMES NB. SITE NAME        TECH. ID.        X/Vx       Y/Vy         Z/Vz         Sigmas                   DATA_START   DATA_END  \n');
    fprintf(fidOffic,'                                            -----------------------m/m/y---------------------  \n');
    fprintf(fidOffic,'----------------------------------------------------------------------------------------------------------------------------------  \n');
end
    %%
idcat=1;


% If the coordinates or velocities were fixed in the global adjustment,
% then write zero as formal error
if globsol.antenna.id_dxyz == 0
   globsol.antenna.sigma_xyz = zeros(nst,3);
end
if globsol.antenna.id_dvxyz == 0
   globsol.antenna.sigma_vxyz = zeros(nst,3);
end


kk=0;
for i=1:nst
    aname=globsol.antenna.refname_coord(i,:);
    stname=cellstr(aname);
    
    dxyz1=globsol.antenna.dxyz(i,:);  % glob mjd
    mxyz1=globsol.antenna.sigma_xyz(i,:);
    
    dvxyz1=globsol.antenna.dvxyz(i,:);
    mvxyz1=globsol.antenna.sigma_vxyz(i,:);

    
    estint=globsol.antenna.epoch_break(i,2:3);
    for j=1:length(refantbr)
        if refantbr(j).name==aname
            ia=j;
            break
        end
    end

    % find apriori intervals included in the estimation intervals
    aprint=refantbr(ia).break_apr;
    
    % intervals beginn at the same time
    id1 = find(estint(1) <= aprint  & aprint  <= estint(2));
    
    if ~isempty(id1)
     epo= refantbr(ia).break_apr(id1);
       
        if isempty(find(estint(1) == aprint));
            id1=[id1(1)-1 id1];
            epo=[ estint(1)  epo];
        end
        if isempty(find(estint(2) == aprint));
            id1=[id1 id1(end)+1 ];
            epo = [epo estint(2)];
        end
    end
    
    if isempty(id1)
        xx = find(estint(1) > aprint );
        id1(1)=xx(end);
        clear xx
        xx = find(estint(2) < aprint );
        id1(2)=xx(1);
        
        epo=estint;
        
    end
    idAPR = id1(1:end-1);
    
    for j=1:length(idAPR)
       X(j)=refantbr(ia).x_apr(idAPR(j)); % mjd a priori catalogue
       Y(j)=refantbr(ia).y_apr(idAPR(j));
       Z(j)=refantbr(ia).z_apr(idAPR(j));
       VX(j)=refantbr(ia).vx_apr(idAPR(j));
       VY(j)=refantbr(ia).vy_apr(idAPR(j));
       VZ(j)=refantbr(ia).vz_apr(idAPR(j));
       
       dxyz(j,:)=dxyz1; % glob mjd
       mxyz(j,:) = mxyz1;
       
       dvxyz(j,:)=dvxyz1;
       mvxyz(j,:) = mvxyz1;
       
       if globsol.antenna.id_dvxyz == 0
           stepoch(j) = refantbr(ia).epoch(idAPR(j));
       else
           stepoch(j) = globsol.antenna.epoch_break(1,1);
       end
       
       mjd_apriori(j) = refantbr(ia).epoch(idAPR(j));
    end

    XYZapr = [X; Y ;Z]';
    Vxyzapr = [VX; VY; VZ]';
    
    XYZapr2mjdglob=[];
    for j=1:length(idAPR)
       XYZapr2mjdglob(j,:) = XYZapr(j,:) +  Vxyzapr(j,:)*(stepoch(j) - mjd_apriori(j))./365.25; %(mjd_glob - mjd_apriori)
    end
    
     Xnew = XYZapr2mjdglob + [dxyz]./100; %mjd_glob
     mXnew = mxyz./100;
     
     
     Vnew = [VX; VY; VZ]'+[dvxyz]./100;
     mVnew = mvxyz./100;
     
     
     for k=1:length(idAPR)
     kk=kk+1;     
         aname_save(kk,:) =  aname;
         Xnew_save(kk,:) =  Xnew(k,:);
         Vnew_save(kk,:) =  Vnew(k,:);
         stepoch_save(kk)= stepoch(k);
         epo_save(kk,:)= [epo(k) epo(k+1)];
     end
     
%         fprintf(fidVieVS,'%c%c%c%c%c%c%c%c     %14.4f  %14.4f  %14.4f      %9.4f   %9.4f   %9.4f      %5.0f   %5.0f   %5.0f \n', ...
%             aname,Xnew(k,:),Vnew(k,:), stepoch(k), epo(k), epo(k+1));
%      end
  
     %%
     
     
    % get index of current station in itrf struct
    aname_ = aname;
    idblank=find(aname_(1:max(find(aname_(:)~=' ')))==' '); %hana Jun11
    aname_(idblank)='_';

    antennaIdx = strcmp({trf.name},aname_);
    if max(antennaIdx)==1
        if sum(antennaIdx)>1
            fprintf('More than one antenna entries match the antenna name')
        end
        antennaIdx=find(antennaIdx);

        %SiteCode
        curSiteCode=str2double(trf(antennaIdx).CDP);
        if isnan(curSiteCode) % if siteCode in File = '----'
            antCAT.siteCode=0;
        else
            antCAT.siteCode=curSiteCode;
        end
        % Domes
        antCAT.domes=trf(antennaIdx).domes;
    else
        antCAT.siteCode=9999;     %hana Jun11
        antCAT.domes=cellstr('99999S000');
    end

    idant = find(strcmp(stname, cellstr(globsol.antenna.antennas)) == 1);

    iact=find(globsol.antenna.antactiv(idant,:)==1);
    antCAT.firstObsMjd = min(globsol.antenna.antactiv(end,iact));
    antCAT.lastObsMjd = max(globsol.antenna.antactiv(end,iact));

    antCAT.firstObsDOY = mjd2yydoysecod(antCAT.firstObsMjd);
    antCAT.lastObsDOY = mjd2yydoysecod(antCAT.lastObsMjd);
        
     for k=1:length(idAPR)
         if epo(k)==0
             firstObsDoy(k,:) = antCAT.firstObsDOY;
             if firstObsDoy(k,3) > 43200
                 firstObsDoy(k,2) = firstObsDoy(k,2)+1; % round the whole days
             end
             firstObsDoy(k,3) = 0;
         else
             firstObsDoy(k,:) = mjd2yydoysecod(epo(k));
         end
         
         if epo(k+1)==99999
             lastObsDoy(k,:) = antCAT.lastObsDOY;
             if lastObsDoy(k,3) > 43200
                 lastObsDoy(k,2) = lastObsDoy(k,2)+1; % round the whole days
             end
             lastObsDoy(k,3) = 0;
         else
             lastObsDoy(k,:) = mjd2yydoysecod(epo(k+1));
         end
        
    end       
    timeYrStrfirst=num2str(firstObsDoy(:,1));
    timeYrStrlast=num2str(lastObsDoy(:,1));
    
    for k=1:length(idAPR)
        catal(idcat).domes = num2str(antCAT.domes);
        catal(idcat).aname = aname_;
        catal(idcat).techn= 'VLBI';
        catal(idcat).siteCode=num2str(antCAT.siteCode);
        catal(idcat).Xnew = Xnew(k,:);
        catal(idcat).mXnew = mXnew(k,:);
        catal(idcat).Vnew = Vnew(k,:);
        catal(idcat).mVnew = mVnew(k,:);
        catal(idcat).timeYrStrfirst =timeYrStrfirst(k,3:4);
        catal(idcat).firstObsDoy = firstObsDoy(k,2:3);
        catal(idcat).timeYrStrlast =timeYrStrlast(k,3:4);
        catal(idcat).lastObsDoy = lastObsDoy(k,2:3);
        
        idcat=idcat+1;
    end
     
     %mVXnew
    clear ic iv iapr i1 inew ia ix iv X Y Z x y z VX VY VZ vx vy vz Xnew Vnew dxyz dvxyz antCAT firstObsDoy lastObsDoy stepoch
    
end

if globsol.antenna.id_dvxyz == 1
    % sort alphabetical
    for i=1:length(catal)
        cataname(i,:) = catal(i).aname;
    end
    [a,idsort]=sortrows(cataname);



     writeFormat1='%9s %8s         %4s %4s %13.4f %13.4f %13.4f   %6.4f %6.4f %6.4f   %2s:%03.0f:%05.0f %2s:%03.0f:%05.0f \n';
     writeFormat2='                                           %7.4f       %7.4f       %7.4f   %6.4f %6.4f %6.4f\n';
     for i=1:length(idsort)
        fprintf(fidOffic,writeFormat1, catal(idsort(i)).domes, catal(idsort(i)).aname, catal(idsort(i)).techn, catal(idsort(i)).siteCode,...
            catal(idsort(i)).Xnew, catal(idsort(i)).mXnew, ...
            catal(idsort(i)).timeYrStrfirst, catal(idsort(i)).firstObsDoy(1), catal(idsort(i)).firstObsDoy(2), ...
            catal(idsort(i)).timeYrStrlast, catal(idsort(i)).lastObsDoy(1), catal(idsort(i)).lastObsDoy(2) );

        fprintf(fidOffic,writeFormat2, catal(idsort(i)).Vnew, catal(idsort(i)).mVnew);

     end
    fclose(fidOffic);
end



%% write a new TRF in a VieVS format

% sort alphabetical
a=[]; idsort=[];
[a,idsort]=sortrows(aname_save);

anameS = aname_save(idsort,:);
XnewS = Xnew_save(idsort,:);
VnewS = Vnew_save(idsort,:);
stepochS = stepoch_save(idsort);
epoS = epo_save(idsort,:);


fidVieVS=fopen([paths.path_out 'TRF/' paths.out '/trf_' paths.L2 '.txt'],'wt');
fprintf(fidVieVS,'%% Created on %02.0f.%02.0f.%04.0f at %02.0f:%02.0f:%02.0f local time\n', curDate(3), curDate(2), curDate(1), curDate(4), curDate(5), curDate(6));
fprintf(fidVieVS,'%% were the station coordinates estimated? (0/1) %1.0f \n',globsol.antenna.id_dxyz);
fprintf(fidVieVS,'%% were the station velocities estimated?  (0/1) %1.0f \n\n\n',globsol.antenna.id_dvxyz);
fprintf(fidVieVS,'%% station          x [m]           y [m]         z [m]            vx [m/y]    vy [m/y]    vz [m/y]      epoch   start    end \n\n');

for k=1:size(anameS,1)
    fprintf(fidVieVS,'%c%c%c%c%c%c%c%c     %14.4f  %14.4f  %14.4f      %9.4f   %9.4f   %9.4f      %5.0f   %5.0f   %5.0f \n', ...
            anameS(k,:),XnewS(k,:),VnewS(k,:), stepochS(k), epoS(k,1), epoS(k,2));
end

fclose(fidVieVS);



