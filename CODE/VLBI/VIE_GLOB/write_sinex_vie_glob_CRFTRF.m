% ************************************************************************
%   Description:
% This function writes SINEX output for the global solution from VIE_GLOB
% based on function write_sinex_vievs.m in VIE_LSM
%
%   Input:										
%       globsol             structure with estimates
%       paths               paths to directories
%       trf                 info from superstation file
%                           
%   Output:                
%      sinex file with N and b, and sinex file with covariance matrix
%
%   External calls: 	
%      mjd2yydoysecod               					    											
%       
%   Coded for VieVS: 
%   23 Oct 2012 by Hana Krásná
%
%   Revision:
%   24 Oct 2012 by Hana Krásná: changes according to change of the units of
%               RA in globsol
%   22 Aug 2013 by Hana Krásná: change in the header
%   18 Dec 2013 by Hana Krásná: ICRF designation for sources added
%   10 Apr 2017 by David Mayer: Added source epochs block
%%
%
% This version ONLY FOR TRF (coord + vel) + CRF (coord) SOLUTION !!!


function write_sinex_vie_glob_CRFTRF(globsol,paths,trf,crf)

    outsnx.soupos=1;
    outsnx.antpos=1;
    outsnx.antvel=1;

    CreateSinex_Cov = 1; % do you want to create a SINEX file with covariance matrix?

    if exist(['../OUT/GLOB/SNX/' paths.out])~=7
        mkdir(['../OUT/GLOB/SNX/' paths.out]);
    end

    snxFile=['../OUT/GLOB/SNX/' paths.out '/' paths.L2 '_Nb.SNX'];
    
    
    

    % define space between two blocks
    blockSpace={'*', '* -----------------------------------------------------------------------------', '*'};

    % open file ('a': writing, append, 'w': discharge content if any)
    fid=fopen(snxFile, 'w');
    
   
    % add site code and itrfXYZ to x_.antenna struct
    nan=size(globsol.antenna.antennas,1);
    for k=1:nan
        % get index of current station in itrf struct
        aname=globsol.antenna.antennas(k,:);
        
        idblank=find(aname(1:max(find(aname(:)~=' ')))==' '); %hana Jun11
        aname(idblank)='_';
        
        a = strfind({trf.name},aname);
        antennaIdx=~cellfun('isempty',a);
        if max(antennaIdx)==1
            if sum(antennaIdx)>1
                fprintf('More than one antenna entries match antenna(k).name')
            end
            antennaIdx=find(antennaIdx);
            
            %SiteCode
            curSiteCode=str2double(trf(antennaIdx).CDP);
            if isnan(curSiteCode) % if siteCode in File = '----'
                antenna(k).siteCode=0;
            else
                antenna(k).siteCode=curSiteCode;
            end
            % Domes
            antenna(k).domes=trf(antennaIdx).domes;
            % PointCode
            antenna(k).pointCode='A';
            %Description (if too long: only 22 characters)
            antenna(k).description=trf(antennaIdx).comments(1:22);
        else
            antenna(k).siteCode=9999;     %hana Jun11
            antenna(k).pointCode='A';
            antenna(k).domes=cellstr('99999S000');
            antenna(k).description=antenna(k).name;
        end
        iact=find(globsol.antenna.antactiv(k,:)==1);
        antenna(k).firstObsMjd=min(globsol.antenna.antactiv(nan+1,iact));
        antenna(k).lastObsMjd=max(globsol.antenna.antactiv(nan+1,iact));
    end

    % add 0 as site codes for stations which do not occur in itrf
    %if sum(cellfun('isempty',{x_.antenna.siteCode})) > 0
    %    x_.antenna(cellfun('isempty',{x_.antenna.siteCode})).siteCode=0;
    %end

    %[pointCode{1:size(x_.antenna, 2), 1}]=deal('A');
    %[obsCode{1:size(x_.antenna, 2), 1}]=deal('R');          % R=VLBI
    %[statDescription{1:size(antenna, 2), 1}]=deal('');
    
    numStat=length(antenna);
    numSou=size(globsol.source.refname.IERS,1);
    
    %% write headerline
    %if existsnx==0
    snxFormat=2.10;         % format of sinex
    creatingAgency='VIE';   % agency creating this file
    curDate=clock;          % current date and time
    doy=datenum(curDate-[curDate(1),0,0,curDate(4:end)]);  % day of year of current time
    secod=curDate(4)*60*60+curDate(5)*60+curDate(6);       % seconds of day of current time
    yrStr=num2str(curDate(1));
    dataAgency='VIE';       % agency providing the data
    mjdsall=globsol.antenna.antactiv(end,:);
    dataStartEnd=mjd2yydoysecod([min(mjdsall), max(mjdsall)]);
%    sessionTimeConsumption=max([antenna.lastObsMjd])-min([antenna.firstObsMjd]);    % needed later when writing epochs where parameters are valid 
    dataStartEndYrStr=num2str(dataStartEnd(:,1));
    obsCode='R';            % observation code (R=VLBI)
    numEstPar=length(globsol.x);
    constrCode=num2str(1);  % 0=tight, 1=significant, 2=unconstrained
    solcontents='S C'; %S - stations param, C - celestial ref. frame 

    %                    2.10 TUW   08:  025 :09862  TUW   00: 000  :00000    00:000   :00000               
    fprintf(fid, '%%=SNX %4.2f %3s %2s:%03.0f:%05.0f %3s %02s:%03.0f:%05.0f %02s:%03.0f:%05.0f %1s %05.0f, %1s, %3s\n', snxFormat, creatingAgency, yrStr(3:4), doy, secod, dataAgency, dataStartEndYrStr(5:2:end), dataStartEnd(1,2), dataStartEnd(1,3), dataStartEndYrStr(6:2:end), dataStartEnd(2,2), dataStartEnd(2,3), obsCode, numEstPar, constrCode, solcontents);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});

    %% write 4 FILE/REFERENCE block if tickbox is ticked
    blockName='FILE/REFERENCE';

    infoType={'DESCRIPTION', 'OUTPUT', 'CONTACT', 'SOFTWARE'};
    info={'Vienna University of Technology', 'Global solution - Vienna VLBI Software', 'hana.krasna@tuwien.ac.at <Hana Krásná>', 'VieVS - Vienna VLBI Software'};
    frBlock=[infoType; info];

    %write data
    fprintf(fid, '+%s\n', blockName);
    fprintf(fid, ' %-18s%-60s\n', frBlock{:});
    fprintf(fid, '-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});
    
    %% write 5 FILE/COMMENT block
    blockName='FILE/COMMENT';
    
    numSes = num2str(size(globsol.sessions,1));
%     numAntDat = num2str(size(globsol.antenna.datum,1));
%     numSouDat = num2str(size(globsol.source.datum,1));
    
    comment={['Input data: ' numSes ' IVS VLBI sessions']};
%     comment={['Input data: ' numSes ' IVS VLBI sessions (1984.0 - 2012.5)'];
%         'Clock offsets, zenit wet delay, tropospheric gradients and EOP';
%         '     were reduced from the session-wise normal equations, i.e. they were estimated implicitly';
%         'TRF (position + velocity) was estimated at epoch J2000.0,';
%         ['     NNT+NNR condition w.r.t. VTRF2008 was applied on ' numAntDat ' stations'];
%         'CRF (position) was estimated,';
%         ['     NNR condition w.r.t. ICRF2 was applied on ' numSouDat ' sources (defining ICRF2)']
%         };
    
    %write data
    fprintf(fid, '+%s\n', blockName);
    fprintf(fid, ' %-79s\n', comment{:});
    fprintf(fid, '-%s\n', blockName);

    
    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});
    
    %% write 11. SOURCE/ID block if tickbox is ticked
    blockName='SOURCE/ID';
    
    % write blockstart line to file
    fprintf(fid, '+%s\n', blockName);
    
    % write description line
    fprintf(fid, '%s\n', '*Code IERS_nam IVS_nam  ICRF_designator');
    
    refnameIERS = globsol.source.refname.IERS;
    refnameIVS = globsol.source.refname.IVS;  
    for i=1:length(crf)
        supsouIERS(i,:)=crf(i).IERSname;
    end

    % for all sources
    for k=1:numSou
        ICRFdes=[''];
        idcrf=(find(strcmp(cellstr(supsouIERS),cellstr(refnameIERS(k,:)))==1));
        if ~isempty(idcrf)
            ICRFdes=crf(idcrf).designation(6:end);
        else
            ICRFdes='J               ';
        end
        fprintf(fid, ' %04s %-8s %-8s %-16s\n', num2str(k), refnameIERS(k,:), refnameIVS(k,:),ICRFdes);
        sources(k).siteCode=k;
    end
    
    % write blockend line to file
    fprintf(fid, '-%s\n', blockName);
    
    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});

    %% write 12. SITE/ID block if tickbox is ticked
    blockName='SITE/ID';
    writeFormat=' %4s %2s %9s %1s %-22s %3.0f %02.0f %04.1f %3.0f %02.0f %04.1f %7.1f\n';

    % get lat and lon from extimated x,y,z
    for i = 1:length(globsol.antenna.apr_pos)
        X(i)=globsol.antenna.apr_pos(i).x_apr(1);
        Y(i)=globsol.antenna.apr_pos(i).y_apr(1);
        Z(i)=globsol.antenna.apr_pos(i).z_apr(1);
        [approxLat,approxLon,approxH]=xyz2ell([X', Y', Z']);
    end
    % rad -> °
    approxLat=approxLat*180/pi;
    approxLon=approxLon*180/pi;

    % [-180, 180] -> [0, 360]
    approxLon(approxLon<0)=approxLon(approxLon<0)+360;

    % ° -> dms
    dmsLat=degrees2dms(approxLat);
    dmsLon=degrees2dms(approxLon);

    % write blockstart line to file
    fprintf(fid, '+%s\n', blockName);

    % write description line
    fprintf(fid, '%s\n', '*CODE PT DOMES____ T STATION_DESCRIPTION___ APPROX_LON_ APPROX_LAT_ APP_H__');
    % write data
    for k=1:length(globsol.antenna.apr_pos)
        %                                  -4s                     %-2s                     %9s             %1s          
        fprintf(fid, writeFormat, num2str(antenna(k).siteCode), antenna(k).pointCode, antenna(k).domes(:), obsCode,...
                antenna(k).description, dmsLon(k,1), dmsLon(k,2), dmsLon(k,3), dmsLat(k,1), dmsLat(k,2), dmsLat(k,3), approxH(k));
              %-22s               %3.0f         %02.0f      %04.1f          %3.0f        %02.0f         %04.1f    7.1f\n';
        %fprintf(fid, '\n');
    end

    % write blockend line to file
    fprintf(fid, '-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});

%     %% write 17. SITE/ECCENTRICITY block if tickbox is ticked
%     blockName='SITE/ECCENTRICITY';
%     writeFormat=' %4.0f %2s %4.0f %1s %2s:%03.0f:%05.0f %2s:%03.0f:%05.0f %3s %8.4f %8.4f %8.4f\n';
%     % write blockstart line to file
%     fprintf(fid, '+%s\n', blockName);
% 
%     % write info line
%     fprintf(fid, '*Code PT SOLN T Data_Start__ Data_End____ typ Up/X____ North/Y_ East/Z__\n');
%     
%     % get number of scans
%     %numScan=size(scan,2);
%     
%     % for all stations
%     for k=1:numStat
%         %curPtCode='A';
%         soln=1; % solution number
%         
%         time=mjd2yydoysecod([antenna(k).firstObsMjd, antenna(k).lastObsMjd, (antenna(k).firstObsMjd+antenna(k).lastObsMjd)/2]);
%         timeYrStr=num2str(time(:,1));
%         
%         % write data
%         fprintf(fid, writeFormat, antenna(k).siteCode, antenna(k).pointCode, soln, obsCode, timeYrStr(7:3:end), time(1,2), time(1,3), timeYrStr(8:3:end), time(2,2), time(2,3), antenna(k).ecctype, antenna(k).c_ecc(1), antenna(k).c_ecc(2), antenna(k).c_ecc(3));
%     end
% 
%     % write blockend line to file
%     fprintf(fid, '-%s\n', blockName);
% 
%     % write headerspace
%     fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});

    %% write 20. SOLUTION/EPOCHS block if tickbox is ticked
    blockName='SOLUTION/EPOCHS';
    writeFormat=' %04s %2s %4.0f %1s %2s:%03.0f:%05.0f %2s:%03.0f:%05.0f %2s:%03.0f:%05.0f\n';

    % write blockstart line to file
    fprintf(fid, '+%s\n', blockName);

    % write info line
    fprintf(fid, '*Code PT SOLN T Data_start__ Data_end____ Mean_epoch__\n');

    % write data
    
    
    % sources   
    souactiv = globsol.source.souactiv;
    grefname = globsol.source.refname;
    
    solID=1;
    pointCodeSou='--';
    for k=1:numSou
        idObsSou=find(souactiv(k,:)>0);
        soumjdall = souactiv(end,idObsSou);
        firstO = min(soumjdall);
        lastO = max(soumjdall);
        
        % get time of first, last and mean mjd of observation
        time=mjd2yydoysecod([firstO, lastO, (firstO+lastO)/2]);
        % year: float 2008 -> string 08 (again for first, last and median scan)
        timeYrStr=num2str(time(:,1));
         
        fprintf(fid, writeFormat, num2str(sources(k).siteCode), pointCodeSou, solID, obsCode, timeYrStr(7:3:end), time(1,2), time(1,3), timeYrStr(8:3:end), time(2,2), time(2,3), timeYrStr(9:3:end), time(3,2), time(3,3));

        clear time firstO lastO idObsSou soumjdall
    end

    
    
    
    % stations
  
    % for each station: one line (since we have only one combination of point
    % code: A / solution ID: 1 / ObservationCode: R=VLBI)
    
    numStatBreak = size(globsol.antenna.refname_coord,1);
    
    for k=1:numStatBreak
        aname = globsol.antenna.refname_coord(k,:);
        if globsol.antenna.epoch_break(k,2)==0
            j=[];
            for i=1:size(globsol.antenna.antennas,1)
                if strcmp(cellstr(aname),cellstr(globsol.antenna.antennas(i,:)))
                    j=i;
                end
            end
            firstO = antenna(j).firstObsMjd;
            solID=1;
        else
            firstO = globsol.antenna.epoch_break(k,2);
            solID=solID+1; %needs to be udpate if more than 1 break occures
        end
        if globsol.antenna.epoch_break(k,3)==99999
            j=[];
            for i=1:size(globsol.antenna.antennas,1)
                if strcmp(cellstr(aname),cellstr(globsol.antenna.antennas(i,:)))
                    j=i;
                end
            end
            lastO = antenna(j).lastObsMjd;
        else
            lastO = globsol.antenna.epoch_break(k,3);
        end
            
        
        % get time of first, last and mean mjd of observation
        time=mjd2yydoysecod([firstO, lastO, (firstO+lastO)/2]);

        % year: float 2008 -> string 08 (again for first, last and median scan)
        timeYrStr=num2str(time(:,1));
        
        %obsCode='R'; % already defined before

        fprintf(fid, writeFormat, num2str(antenna(j).siteCode), antenna(j).pointCode, solID, obsCode, timeYrStr(7:3:end), time(1,2), time(1,3), timeYrStr(8:3:end), time(2,2), time(2,3), timeYrStr(9:3:end), time(3,2), time(3,3));
        clear time
    end

    % write blockend line to file
    fprintf(fid, '-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});
    %% write SOURCE/EPOCHS block ##########################
    if outsnx.soupos
        blockName='SOURCE/EPOCHS';
        writeFormat=' %04.0f       %1.0f %1s %2s:%03.0f:%05.0f %2s:%03.0f:%05.0f %2s:%03.0f:%05.0f %10.0f %10.0f\n';
        % write blockstart line to file
        fprintf(fid, '+%s\n', blockName);
        % write info line
        fprintf(fid, '*Code     SID T Data_start__ Data_end____ Mean_epoch__ Num_of_Ses Num_of_Obs\n');
        
        % sources
        souactiv = globsol.source.souactiv;
        
        solID=1;
        pointCodeSou='--';
        for k=1:numSou
            idObsSou=find(souactiv(k,:)>0);
            soumjdall = souactiv(end,idObsSou);
            
            source_numSess = length(idObsSou);
            source_numobs = sum(souactiv(k,idObsSou));
            
            firstO = min(soumjdall);
            lastO = max(soumjdall);
            
            % get time of first, last and mean mjd of observation
            time=mjd2yydoysecod([firstO, lastO, (firstO+lastO)/2]);
            % year: float 2008 -> string 08 (again for first, last and median scan)
            timeYrStr=num2str(time(:,1));
            fprintf(fid, writeFormat, sources(k).siteCode, solID, obsCode, timeYrStr(7:3:end), time(1,2), time(1,3), timeYrStr(8:3:end), time(2,2), time(2,3), timeYrStr(9:3:end), time(3,2), time(3,3),source_numSess,source_numobs);

            %fprintf(fid, writeFormat, num2str(sources(k).siteCode), pointCodeSou, solID, obsCode, timeYrStr(7:3:end), time(1,2), time(1,3), timeYrStr(8:3:end), time(2,2), time(2,3), timeYrStr(9:3:end), time(3,2), time(3,3),source_numSess,source_numobs);
            
            clear time firstO lastO idObsSou soumjdall
        end
        % write blockend line to file
        fprintf(fid, '-%s\n', blockName);
        
        % write headerspace
        fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});
    end
    
    %% write 22 SOLUTION/STATISTICS block
    % define name of block and output format
    blockName='SOLUTION/STATISTICS';
    writeFormat=' %-30s %22.15f\n';
    
   
    %if estimation was carried out
    infoType={'NUMBER OF OBSERVATIONS', 'NUMBER OF UNKNOWNS', 'SQUARE SUM OF RESIDUALS (VTPV)', 'VARIANCE FACTOR' }; % hana
    info={globsol.numObs, globsol.numEst, globsol.vTPv, globsol.varfac }; % hana

    frBlock=[infoType; info];

    %write data
    fprintf(fid, '+%s\n', blockName);
    fprintf(fid, writeFormat, frBlock{:});
    fprintf(fid, '-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});
        
    %% write 23 SOLUTION/ESTIMATES block if tickbox is ticked
    
    aprDate=mjd2yydoysecod(globsol.antenna.epoch_break(1,1));
    aprDateYrStr=num2str(aprDate(1));
    
    % if option was chosen in gui
        
        % define name of block and output format
        blockName='SOLUTION/ESTIMATES';
        writeFormat=' %5.0f %-6s %04s %2s %4s %2s:%03.0f:%05.0f %-4s %1.0f %21s %11s\n';
        %writeFormat=' %5.0f %-6s %4s %2s %4s %2s:%03.0f:%05.0f %-4s %1.0f %21.13d %11.4d\n';
        if ispc, formatEstVal='%22.14e'; formatStDev='%12.5e'; else
                 formatEstVal='%21.14e'; formatStDev='%11.5e';
        end
        constrEstCoord=1;
        curIndex=1;

        % write blockstart line to file
        fprintf(fid, '+%s\n', blockName);

        % write info line
        fprintf(fid, '*Index Type__ Code Pt Soln Ref_Epoch___ Unit S Estimated_Value______ Std_Dev____\n');

        % write data  

         if outsnx.soupos==1
            glsou = globsol.source;
            for k=1:numSou
                % Solution ID   
                soln=num2str(1);

                % calculate total estimated values
                sources(k).totRa = glsou.apriori_rade(k,1) + glsou.drade(k,1)/1000/3600/12*pi; % ms --> rad
                sources(k).totDe = glsou.apriori_rade(k,2) + glsou.drade(k,2)/1000/3600/180*pi; % mas --> rad
                sources(k).totRaStd = glsou.sigma_rade(k,1)/1000/3600/12*pi; 
                sources(k).totDeStd = glsou.sigma_rade(k,2)/1000/3600/180*pi; 
                constrEstSou=1;

                % preparation of values for format e+02 instead of e+002
                totRastring=sprintf(formatEstVal, sources(k).totRa);
                if ispc, totRastring = strrep(totRastring, 'e+0', 'e+'); totRastring = strrep(totRastring, 'e-0', 'e-');   end
                totDestring=sprintf(formatEstVal, sources(k).totDe);
                if ispc, totDestring = strrep(totDestring, 'e+0', 'e+'); totDestring = strrep(totDestring, 'e-0', 'e-'); end

                totRastdString=sprintf(formatStDev, sources(k).totRaStd);
                if ispc, totRastdString = strrep(totRastdString, 'e+0', 'e+'); totRastdString = strrep(totRastdString, 'e-0', 'e-');   end
                totDestdString=sprintf(formatStDev,sources(k).totDeStd);
                if ispc, totDestdString = strrep(totDestdString, 'e+0', 'e+'); totDestdString = strrep(totDestdString, 'e-0', 'e-');  end

                fprintf(fid, writeFormat, curIndex,  'RS_RA', num2str(sources(k).siteCode), '--', soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'rad', constrEstSou, totRastring, totRastdString);
                fprintf(fid, writeFormat, curIndex+1,'RS_DE', num2str(sources(k).siteCode), '--', soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'rad', constrEstSou, totDestring, totDestdString);

                curIndex=curIndex+2;
            end
        end

        
        
        if outsnx.antpos==1
            % write x,y,z coordinates
            glbapr=globsol.antenna.apr_pos;
            u=1;
            for k=1:numStat

                intid=find(globsol.antenna.apr_pos(k).interv>0);
                for j=1:length(intid)
                fi = [];
                % Solution ID   
                soln=num2str(j);

                breakap = glbapr(k).break(intid(j));
                fi=find(breakap == glbapr(k).break_apr);
                if ~isempty(fi)
                    idc = fi;
                else
                    fi = find(breakap > glbapr(k).break_apr);
                    idc=fi(end);
                end
                
                dxyz=globsol.antenna.dxyz(u,:);
                
                dvxyz=[0 0 0];
                if outsnx.antvel==1
                    dvxyz=globsol.antenna.dvxyz(u,:);
                end
                
                % calculate apriori values
                %                  coordx   +     vx      * (        mean scan mjd                          - itrf epoch mjd )/ diff(mjd)->years                                
                antenna(k).aprX=glbapr(k).x_apr(idc);
                antenna(k).aprY=glbapr(k).y_apr(idc);
                antenna(k).aprZ=glbapr(k).z_apr(idc);

                % velocity
                antenna(k).aprVX=glbapr(k).vx_apr(idc);
                antenna(k).aprVY=glbapr(k).vy_apr(idc);
                antenna(k).aprVZ=glbapr(k).vz_apr(idc);
                
                % calculate total estimated values
                antenna(k).totX=antenna(k).aprX+dxyz(1)/100;
                antenna(k).totY=antenna(k).aprY+dxyz(2)/100;
                antenna(k).totZ=antenna(k).aprZ+dxyz(3)/100;
                
                antenna(k).totVX=antenna(k).aprVX+dvxyz(1)/100;
                antenna(k).totVY=antenna(k).aprVY+dvxyz(2)/100;
                antenna(k).totVZ=antenna(k).aprVZ+dvxyz(3)/100;

                
                % preparation of values for format e+02 instead of e+002
                totXstring=sprintf(formatEstVal, antenna(k).totX);
                if ispc, totXstring = strrep(totXstring, 'e+0', 'e+'); totXstring = strrep(totXstring, 'e-0', 'e-');   end
                totYstring=sprintf(formatEstVal, antenna(k).totY);
                if ispc, totYstring = strrep(totYstring, 'e+0', 'e+'); totYstring = strrep(totYstring, 'e-0', 'e-'); end
                totZstring=sprintf(formatEstVal, antenna(k).totZ);
                if ispc, totZstring = strrep(totZstring, 'e+0', 'e+'); totZstring = strrep(totZstring, 'e-0', 'e-'); end

                totVXstring=sprintf(formatEstVal, antenna(k).totVX);
                if ispc, totVXstring = strrep(totVXstring, 'e+0', 'e+'); totVXstring = strrep(totVXstring, 'e-0', 'e-');   end
                totVYstring=sprintf(formatEstVal, antenna(k).totVY);
                if ispc, totVYstring = strrep(totVYstring, 'e+0', 'e+'); totVYstring = strrep(totVYstring, 'e-0', 'e-'); end
                totVZstring=sprintf(formatEstVal, antenna(k).totVZ);
                if ispc, totVZstring = strrep(totVZstring, 'e+0', 'e+'); totVZstring = strrep(totVZstring, 'e-0', 'e-'); end

                
                
                totXstdString=sprintf(formatStDev, globsol.antenna.sigma_xyz(u,1)/100);
                if ispc, totXstdString = strrep(totXstdString, 'e+0', 'e+'); totXstdString = strrep(totXstdString, 'e-0', 'e-');   end
                totYstdString=sprintf(formatStDev, globsol.antenna.sigma_xyz(u,2)/100);
                if ispc, totYstdString = strrep(totYstdString, 'e+0', 'e+'); totYstdString = strrep(totYstdString, 'e-0', 'e-');  end
                totZstdString=sprintf(formatStDev, globsol.antenna.sigma_xyz(u,3)/100);
                if ispc, totZstdString = strrep(totZstdString, 'e+0', 'e+'); totZstdString = strrep(totZstdString, 'e-0', 'e-');  end
                
                fprintf(fid, writeFormat, curIndex,  'STAX', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrEstCoord, totXstring, totXstdString);
                fprintf(fid, writeFormat, curIndex+1,'STAY', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrEstCoord, totYstring, totYstdString);
                fprintf(fid, writeFormat, curIndex+2,'STAZ', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrEstCoord, totZstring, totZstdString);
               
                addIndAnt=3; % only pos
                if outsnx.antvel==1
                    totVXstdString=sprintf(formatStDev, globsol.antenna.sigma_vxyz(u,1)/100);
                    if ispc, totVXstdString = strrep(totVXstdString, 'e+0', 'e+'); totVXstdString = strrep(totVXstdString, 'e-0', 'e-');   end
                    totVYstdString=sprintf(formatStDev, globsol.antenna.sigma_vxyz(u,2)/100);
                    if ispc, totVYstdString = strrep(totVYstdString, 'e+0', 'e+'); totVYstdString = strrep(totVYstdString, 'e-0', 'e-');  end
                    totVZstdString=sprintf(formatStDev, globsol.antenna.sigma_vxyz(u,3)/100);
                    if ispc, totVZstdString = strrep(totVZstdString, 'e+0', 'e+'); totVZstdString = strrep(totVZstdString, 'e-0', 'e-');  end


                    fprintf(fid, writeFormat, curIndex+3,'VELX', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm/y', constrEstCoord, totVXstring, totVXstdString);
                    fprintf(fid, writeFormat, curIndex+4,'VELY', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm/y', constrEstCoord, totVYstring, totVYstdString);
                    fprintf(fid, writeFormat, curIndex+5,'VELZ', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm/y', constrEstCoord, totVZstring, totVZstdString);
                
                    addIndAnt = 6; % pos + vel
                end
                
             curIndex=curIndex+addIndAnt;   
             u=u+1;
                end
            end
        end
        
        % write blockend line to file
        fprintf(fid, '-%s\n', blockName);

        % write headerspace
        fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});
    

    %% write 24. SOLUTION/APRIORI block if tickbox is ticked
    blockName='SOLUTION/APRIORI';
    writeFormat=' %5.0f %-6s %04s %2s %4s %2s:%03.0f:%05.0f %-4s %1.0f %21s %11s\n';
    if ispc, formatAprVal='%22.14e'; formatStDev='%12.5e'; else
             formatAprVal='%21.14e'; formatStDev='%11.5e';
    end

    curIndex=1;
    constrApr=2;
    
    % write blockstart line to file
    fprintf(fid, '+%s\n', blockName);

    % write info line
    fprintf(fid, '*Index Type__ CODE PT Soln Ref_epoch___ Unit S Apriori_value________ Std_Dev____\n');

    % write data
    
    % ra, de sources
    if outsnx.soupos==1    
        for k=1:numSou
            % solution ID
            soln=num2str(1);


            sources(k).aprra_sigma=0;
            sources(k).aprde_sigma=0;

            % get format e+02
            aprra=sprintf(formatAprVal, glsou.apriori_rade(k,1));
            if ispc, aprra = strrep(aprra, 'e+0', 'e+'); aprra = strrep(aprra, 'e-0', 'e-');   end
            aprde=sprintf(formatAprVal, glsou.apriori_rade(k,2));
            if ispc, aprde = strrep(aprde, 'e+0', 'e+'); aprde = strrep(aprde, 'e-0', 'e-');  end

            aprra_sigma=sprintf(formatStDev, sources(k).aprra_sigma);
            if ispc, aprra_sigma = strrep(aprra_sigma, 'e+0', 'e+'); aprra_sigma = strrep(aprra_sigma, 'e-0', 'e-');  end
            aprde_sigma=sprintf(formatStDev, sources(k).aprde_sigma);
            if ispc, aprde_sigma = strrep(aprde_sigma, 'e+0', 'e+'); aprde_sigma = strrep(aprde_sigma, 'e-0', 'e-');   end

            % write one line for each RA De
            fprintf(fid, writeFormat, curIndex,  'RS_RA',num2str(sources(k).siteCode), '--', soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'rad', constrApr, aprra, aprra_sigma);
            fprintf(fid, writeFormat, curIndex+1,'RS_DE',num2str(sources(k).siteCode), '--', soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'rad', constrApr, aprde, aprde_sigma);
            curIndex=curIndex+2;
        end
    end
    
    if outsnx.antpos==1
        % x,y,z coordinates
        glbapr=globsol.antenna.apr_pos;
        u=1;
        for k=1:numStat
            intid=find(globsol.antenna.apr_pos(k).interv>0);
            for j=1:length(intid)
                fi = [];
                % Solution ID   
                soln=num2str(j);

                breakap = glbapr(k).break(intid(j));
                fi=find(breakap == glbapr(k).break_apr);
                if ~isempty(fi)
                    idc = fi;
                else
                    fi = find(breakap > glbapr(k).break_apr);
                    idc=fi(end);
                end
                
                
                % calculate apriori values
                antenna(k).aprX=glbapr(k).x_apr(idc);
                antenna(k).aprY=glbapr(k).y_apr(idc);
                antenna(k).aprZ=glbapr(k).z_apr(idc);
                % velocity
                antenna(k).aprVX=glbapr(k).vx_apr(idc);
                antenna(k).aprVY=glbapr(k).vy_apr(idc);
                antenna(k).aprVZ=glbapr(k).vz_apr(idc);

            
                antenna(k).aprX_sigma=0;
                antenna(k).aprY_sigma=0;
                antenna(k).aprZ_sigma=0;
                antenna(k).aprVX_sigma=0;
                antenna(k).aprVY_sigma=0;
                antenna(k).aprVZ_sigma=0;

                % get format e+02
                aprX=sprintf(formatAprVal, antenna(k).aprX);
                if ispc, aprX = strrep(aprX, 'e+0', 'e+'); aprX = strrep(aprX, 'e-0', 'e-');  end
                aprY=sprintf(formatAprVal, antenna(k).aprY);
                if ispc, aprY = strrep(aprY, 'e+0', 'e+'); aprY = strrep(aprY, 'e-0', 'e-');  end
                aprZ=sprintf(formatAprVal, antenna(k).aprZ);
                if ispc, aprZ = strrep(aprZ, 'e+0', 'e+'); aprZ = strrep(aprZ, 'e-0', 'e-');  end

                if outsnx.antvel==1
                    aprVX=sprintf(formatAprVal, antenna(k).aprVX);
                    if ispc, aprVX = strrep(aprVX, 'e+0', 'e+'); aprVX = strrep(aprVX, 'e-0', 'e-');  end
                    aprVY=sprintf(formatAprVal, antenna(k).aprVY);
                    if ispc, aprVY = strrep(aprVY, 'e+0', 'e+'); aprVY = strrep(aprVY, 'e-0', 'e-');  end
                    aprVZ=sprintf(formatAprVal, antenna(k).aprVZ);
                    if ispc, aprVZ = strrep(aprVZ, 'e+0', 'e+'); aprVZ = strrep(aprVZ, 'e-0', 'e-');  end
                end
                
                aprX_sigma=sprintf(formatStDev, antenna(k).aprX_sigma);
                if ispc, aprX_sigma = strrep(aprX_sigma, 'e+0', 'e+'); aprX_sigma = strrep(aprX_sigma, 'e-0', 'e-');  end
                aprY_sigma=sprintf(formatStDev, antenna(k).aprY_sigma);
                if ispc, aprY_sigma = strrep(aprY_sigma, 'e+0', 'e+'); aprY_sigma = strrep(aprY_sigma, 'e-0', 'e-');   end
                aprZ_sigma=sprintf(formatStDev, antenna(k).aprZ_sigma);
                if ispc, aprZ_sigma = strrep(aprZ_sigma, 'e+0', 'e+'); aprZ_sigma = strrep(aprZ_sigma, 'e-0', 'e-');  end

                if outsnx.antvel==1
                    aprVX_sigma=sprintf(formatStDev, antenna(k).aprVX_sigma);
                    if ispc, aprVX_sigma = strrep(aprVX_sigma, 'e+0', 'e+'); aprVX_sigma = strrep(aprVX_sigma, 'e-0', 'e-');  end
                    aprVY_sigma=sprintf(formatStDev, antenna(k).aprVY_sigma);
                    if ispc, aprVY_sigma = strrep(aprVY_sigma, 'e+0', 'e+'); aprVY_sigma = strrep(aprVY_sigma, 'e-0', 'e-');   end
                    aprVZ_sigma=sprintf(formatStDev, antenna(k).aprVZ_sigma);
                    if ispc, aprVZ_sigma = strrep(aprVZ_sigma, 'e+0', 'e+'); aprVZ_sigma = strrep(aprVZ_sigma, 'e-0', 'e-');  end
                end

                % write one line for each x,y,z
                fprintf(fid, writeFormat, curIndex,  'STAX', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrApr, aprX, aprX_sigma);
                fprintf(fid, writeFormat, curIndex+1,'STAY', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrApr, aprY, aprY_sigma);
                fprintf(fid, writeFormat, curIndex+2,'STAZ', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrApr, aprZ, aprZ_sigma);
                addIndAnt=3;
                
                if outsnx.antvel==1
                    fprintf(fid, writeFormat, curIndex+3,'VELX', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm/y', constrApr, aprVX, aprVX_sigma);
                    fprintf(fid, writeFormat, curIndex+4,'VELY', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm/y', constrApr, aprVY, aprVY_sigma);
                    fprintf(fid, writeFormat, curIndex+5,'VELZ', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm/y', constrApr, aprVZ, aprVZ_sigma);
                    addIndAnt=6;
                end

                curIndex=curIndex+addIndAnt;
            end
        end
    end
    
    

    % write blockend line to file
    fprintf(fid, '-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});

    %% change units of the N matrix, b vector and the covariance matrix
    N_sinex = globsol.Nfree;
    b_sinex = globsol.bfree;
    Q_sinex = globsol.Q(1:curIndex-1,1:curIndex-1);
    
    c_xyz = [globsol.antenna.dxyzcol(:,1); globsol.antenna.dxyzcol(:,2);  globsol.antenna.dxyzcol(:,3)];
    c_vxyz = [globsol.antenna.dvxyzcol(:,1); globsol.antenna.dvxyzcol(:,2);  globsol.antenna.dvxyzcol(:,3)];
    c_sou = [globsol.source.dradecol(:,1); globsol.source.dradecol(:,2)];
    
    N_sinex(c_xyz,c_xyz) = N_sinex(c_xyz,c_xyz).*10000; %1/cm^2 --> 1/m^2
    N_sinex(c_xyz,c_vxyz) = N_sinex(c_xyz,c_vxyz).*(100*100*100); %100y/cm^2 --> y/m^2
    N_sinex(c_xyz,c_sou) = N_sinex(c_xyz,c_sou).*(100*1000*3600*180/pi); %1/mas.cm --> 1/rad.m
    N_sinex(c_vxyz,c_xyz) = N_sinex(c_vxyz,c_xyz).*(100*100*100); %100y/cm^2 --> y/m^2 
    N_sinex(c_vxyz,c_vxyz) = N_sinex(c_vxyz,c_vxyz).*(10000*10000); %10000y^2/cm^2 --> y^2/m^2
    N_sinex(c_vxyz,c_sou) = N_sinex(c_vxyz,c_sou).*(100*1000*3600*180/pi*100); %100y/mas.cm --> y/rad.m
    N_sinex(c_sou,c_xyz) = N_sinex(c_sou,c_xyz).*(100*1000*3600*180/pi); %1/mas.cm --> 1/rad.m
    N_sinex(c_sou,c_vxyz) = N_sinex(c_sou,c_vxyz).*(100*1000*3600*180/pi*100); %100y/mas.cm --> y/rad.m
    N_sinex(c_sou,c_sou) = N_sinex(c_sou,c_sou).*(1000*3600*180/pi)^2; %1/mas^2 --> 1/rad^2
    
    b_sinex(c_xyz)=b_sinex(c_xyz).*100; %1/cm --> 1/m
    b_sinex(c_vxyz)=b_sinex(c_vxyz).*(100*100); %100y/cm --> y/m
    b_sinex(c_sou)=b_sinex(c_sou).*(1000*3600*180/pi); %1/mas --> 1/rad
    
    Q_sinex(c_xyz,c_xyz) = Q_sinex(c_xyz,c_xyz)./10000; %cm^2 --> m^2
    Q_sinex(c_xyz,c_vxyz) = Q_sinex(c_xyz,c_vxyz)./(100*100*100); %cm^2/100y --> m^2/y
    Q_sinex(c_xyz,c_sou) = Q_sinex(c_xyz,c_sou)./(100*1000*3600*180/pi); %mas.cm --> rad.m
    Q_sinex(c_vxyz,c_xyz) = Q_sinex(c_vxyz,c_xyz)./(100*100*100); %cm^2/100y --> m^2/y
    Q_sinex(c_vxyz,c_vxyz) = Q_sinex(c_vxyz,c_vxyz)./(10000*10000); %cm^2/10000y^2 --> m^2/y^2
    Q_sinex(c_vxyz,c_sou) = Q_sinex(c_vxyz,c_sou)./(100*1000*3600*180/pi*100); %mas.cm/100y --> rad.m/y
    Q_sinex(c_sou,c_xyz) = Q_sinex(c_sou,c_xyz)./(100*1000*3600*180/pi); %mas.cm --> rad.m
    Q_sinex(c_sou,c_vxyz) = Q_sinex(c_sou,c_vxyz)./(100*1000*3600*180/pi*100); %mas.cm/100y --> rad.m/y
    Q_sinex(c_sou,c_sou) = Q_sinex(c_sou,c_sou)./(1000*3600*180/pi)^2; %mas^2 --> rad^2

    
    %%
    col_x = globsol.antenna.dxyzcol(:,1);
    col_y = globsol.antenna.dxyzcol(:,2);
    col_z = globsol.antenna.dxyzcol(:,3);
    col_vx = globsol.antenna.dvxyzcol(:,1);
    col_vy = globsol.antenna.dvxyzcol(:,2);
    col_vz = globsol.antenna.dvxyzcol(:,3);
    col_ra = globsol.source.dradecol(:,1);
    col_de = globsol.source.dradecol(:,2);
    
    
    if CreateSinex_Cov ==1    
        copyfile(snxFile,[snxFile(:,1:end-7) '_Cov.SNX' ])
    end
    
    %% write 27. SOLUTION/NORMAL_EQUATION_VECTOR block if tickbox is ticked
    blockName='SOLUTION/NORMAL_EQUATION_VECTOR';
    writeFormat=' %5.0f %-6s %04s %2s %4s %2s:%03.0f:%05.0f %4s %1.0f %21s\n';
    if ispc, formatVectorValue='%22.14e'; else
             formatVectorValue='%21.14e'; 
    end
    constrNevCoord=2;
    % define varible index (one for each estimate)
    curIndex=1;

    % write blockstart line to file
    fprintf(fid, '+%s\n', blockName);

    % write info line
    fprintf(fid, '*Index Type__ Code Pt Soln Ref_Epoch___ Unit S RightHandSideVector_b\n');

    % write b vector to file


    % write b vector for sources
    if outsnx.soupos==1
        constrNevSou=2;
        for k=1:numSou
            % make proper format (e+05 instead of e+005)
            ra_b=sprintf(formatVectorValue, b_sinex(col_ra(k)));
            if ispc, ra_b = strrep(ra_b, 'e+0', 'e+'); ra_b = strrep(ra_b, 'e-0', 'e-');  end
            de_b=sprintf(formatVectorValue, b_sinex(col_de(k)));
            if ispc, de_b = strrep(de_b, 'e+0', 'e+'); de_b = strrep(de_b, 'e-0', 'e-');   end

            fprintf(fid, writeFormat, curIndex,  'RS_RA', num2str(sources(k).siteCode), '--', soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'rad', constrNevSou, ra_b);
            fprintf(fid, writeFormat, curIndex+1,'RS_DE', num2str(sources(k).siteCode), '--', soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'rad', constrNevSou, de_b);
            curIndex=curIndex+2;
        end
    end

    % write b vector for station coordinates
    if outsnx.antpos==1
        u=1;
        for k=1:numStat
            intid=find(globsol.antenna.apr_pos(k).interv>0);
            for j=1:length(intid)
                % Solution ID   
                soln=num2str(j);


                % make proper format (e+05 instead of e+005)
                stax_b=sprintf(formatVectorValue, b_sinex(col_x(u)));
                if ispc, stax_b = strrep(stax_b, 'e+0', 'e+'); stax_b = strrep(stax_b, 'e-0', 'e-');  end
                stay_b=sprintf(formatVectorValue, b_sinex(col_y(u)));
                if ispc, stay_b = strrep(stay_b, 'e+0', 'e+'); stay_b = strrep(stay_b, 'e-0', 'e-');   end
                staz_b=sprintf(formatVectorValue, b_sinex(col_z(u)));
                if ispc, staz_b = strrep(staz_b, 'e+0', 'e+'); staz_b = strrep(staz_b, 'e-0', 'e-');   end

                stavx_b=sprintf(formatVectorValue, b_sinex(col_vx(u)));
                if ispc, stavx_b = strrep(stavx_b, 'e+0', 'e+'); stavx_b = strrep(stavx_b, 'e-0', 'e-');  end
                stavy_b=sprintf(formatVectorValue, b_sinex(col_vy(u)));
                if ispc, stavy_b = strrep(stavy_b, 'e+0', 'e+'); stavy_b = strrep(stavy_b, 'e-0', 'e-');   end
                stavz_b=sprintf(formatVectorValue, b_sinex(col_vz(u)));
                if ispc, stavz_b = strrep(stavz_b, 'e+0', 'e+'); stavz_b = strrep(stavz_b, 'e-0', 'e-');   end


                fprintf(fid, writeFormat, curIndex, 'STAX',   num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrNevCoord, stax_b);
                fprintf(fid, writeFormat, curIndex+1, 'STAY', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrNevCoord, stay_b);
                fprintf(fid, writeFormat, curIndex+2, 'STAZ', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrNevCoord, staz_b);

                fprintf(fid, writeFormat, curIndex+3, 'VELX',   num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm/y', constrNevCoord, stavx_b);
                fprintf(fid, writeFormat, curIndex+4, 'VELY', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm/y', constrNevCoord, stavy_b);
                fprintf(fid, writeFormat, curIndex+5, 'VELZ', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm/y', constrNevCoord, stavz_b);

                curIndex=curIndex+6;
                u=u+1;
            end
        end
    end


    % write blockend line to file
    fprintf(fid, '-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});



    %% write 28. SOLUTION/NORMAL_EQUATION_MATRIX block if tickbox is ticked
    blockName='SOLUTION/NORMAL_EQUATION_MATRIX L';
    %writeFormat=' %5.0f %5.0f %s\n';
    formatNormEquInd=' %5.0f %5.0f';
    if ispc, formatNormEquVal='  %21.13e';  else
             formatNormEquVal='  %20.13e';
    end

    % write blockstart line to file
    fprintf(fid, '+%s\n', blockName);

    % write info line
    fprintf(fid, '*Row__ Col__ Norm_Equ_Matrix_Value Norm_Equ_Matrix_Valu2 Norm_Equ_Matrix_Valu3');


    % get all indices for N Matrix in one array
    tmpMatCoorVel=[globsol.antenna.dxyzcol globsol.antenna.dvxyzcol]';
    %tmpMatVel=globsol.antenna.dvxyzcol';
    tmpMatSou=globsol.source.dradecol';

    % make one vector out of it...
    tmpMat=[tmpMatSou(:); tmpMatCoorVel(:)];

    % ... and delete zeros
    tmpMat(tmpMat==0)=[];

    % reorder the normal equation matrix
    % preallocate
    N=zeros(size(N_sinex, 1), size(N_sinex, 2));


    numStatEst=numStatBreak;
%     if outsnx.antpos==0; numStatEst=0; end % if station coordinates are fixed

    for col=1:size(N,1)

        for sou=1:numSou
            N((sou-1)*2+1, col)=N_sinex(col_ra(sou), tmpMat(col));
            N((sou-1)*2+2, col)=N_sinex(col_de(sou), tmpMat(col));
        end

        for stat=1:numStatEst
            N(numSou*2+(stat-1)*6+1,col)=N_sinex(col_x(stat), tmpMat(col));
            N(numSou*2+(stat-1)*6+2,col)=N_sinex(col_y(stat), tmpMat(col));
            N(numSou*2+(stat-1)*6+3,col)=N_sinex(col_z(stat), tmpMat(col));

            N(numSou*2+(stat-1)*6+4,col)=N_sinex(col_vx(stat), tmpMat(col));
            N(numSou*2+(stat-1)*6+5,col)=N_sinex(col_vy(stat), tmpMat(col));
            N(numSou*2+(stat-1)*6+6,col)=N_sinex(col_vz(stat), tmpMat(col));
        end


    end

    % write normal equation matrix
    for row=1:size(N,1)
        for col=1:row
            if mod(col,3)==1
                fprintf(fid, '\n');
                fprintf(fid, formatNormEquInd, row, col);
            end
            % make proper format
            Nstring=sprintf(formatNormEquVal, N(row,col));
            if ispc, Nstring = strrep(Nstring, 'e+0', 'e+'); Nstring = strrep(Nstring, 'e-0', 'e-');   end
               fprintf(fid, '%21s', Nstring);
        end
    end


    % write blockend line to file
    fprintf(fid, '\n-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});

    %% write endsinex line
    fprintf(fid, '%%ENDSNX');
    fclose(fid);    
    
    
    if CreateSinex_Cov ==1
        fidCov = fopen([snxFile(:,1:end-7) '_Cov.SNX'], 'a');
    
        %% write 25. SOLUTION/MATRIX_ESTIMATE block if tickbox is ticked
        blockName='SOLUTION/MATRIX_ESTIMATE L COVA';
        %writeFormat=' %5.0f %5.0f %s\n';
        formatNormEquInd=' %5.0f %5.0f';
        if ispc, formatNormEquVal='  %21.13e';  else
                 formatNormEquVal='  %20.13e';
        end

        % write blockstart line to file
        fprintf(fidCov, '+%s\n', blockName);

        % write info line
        fprintf(fidCov, '*Row__ Col__  Matrix_Estim_Value    Matrix_Estim_Valu2    Matrix_Estim_Valu3');


        % get all indices for N Matrix in one array
        tmpMatCoorVel=[globsol.antenna.dxyzcol globsol.antenna.dvxyzcol]';
        %tmpMatVel=globsol.antenna.dvxyzcol';
        tmpMatSou=globsol.source.dradecol';

        % make one vector out of it...
        tmpMat=[tmpMatSou(:); tmpMatCoorVel(:)];

        % ... and delete zeros
        tmpMat(tmpMat==0)=[];

        % reorder the normal equation matrix
        % preallocate
        Q=zeros(size(Q_sinex, 1), size(Q_sinex, 2));


        numStatEst=numStatBreak;
    %     if outsnx.antpos==0; numStatEst=0; end % if station coordinates are fixed

        for col=1:size(Q,1)

            for sou=1:numSou
                Q((sou-1)*2+1, col)=Q_sinex(col_ra(sou), tmpMat(col));
                Q((sou-1)*2+2, col)=Q_sinex(col_de(sou), tmpMat(col));
            end

            for stat=1:numStatEst
                Q(numSou*2+(stat-1)*6+1,col)=Q_sinex(col_x(stat), tmpMat(col));
                Q(numSou*2+(stat-1)*6+2,col)=Q_sinex(col_y(stat), tmpMat(col));
                Q(numSou*2+(stat-1)*6+3,col)=Q_sinex(col_z(stat), tmpMat(col));

                Q(numSou*2+(stat-1)*6+4,col)=Q_sinex(col_vx(stat), tmpMat(col));
                Q(numSou*2+(stat-1)*6+5,col)=Q_sinex(col_vy(stat), tmpMat(col));
                Q(numSou*2+(stat-1)*6+6,col)=Q_sinex(col_vz(stat), tmpMat(col));
            end


        end

        % write normal equation matrix
        for row=1:size(Q,1)
            for col=1:row
                if mod(col,3)==1
                    fprintf(fidCov, '\n');
                    fprintf(fidCov, formatNormEquInd, row, col);
                end
                % make proper format
                Nstring=sprintf(formatNormEquVal, Q(row,col));
                if ispc, Nstring = strrep(Nstring, 'e+0', 'e+'); Nstring = strrep(Nstring, 'e-0', 'e-');   end
                   fprintf(fidCov, '%21s', Nstring);
            end
        end


        % write blockend line to file
        fprintf(fidCov, '\n-%s\n', blockName);

        % write headerspace
        fprintf(fidCov, '%s\n%s\n%s\n', blockSpace{:});
        %% write endsinex line
        fprintf(fidCov, '%%ENDSNX');
        fclose(fidCov);
    end

    
    %% delete variables
    clear obsCode
    
    %fprintf('%2.0f', pl)
