% ************************************************************************
%   Description:
%   This function exports simulated group delay observables to ascii files
%   in NGS file format. The value actually written to data card # 2 is
%   (o-c) + c - corr, where (o-c) is simulated. And corr includes cable cal
%   and inosphere correction. This means that after reading the simulated
%   NGS files with read_ngs.m the value in scan.obs.obs will be exactly the
%   simulated (o-c)+c. vie_lsm will then compute the O-C as usual, which
%   means that exactly the simulated (o-c) will enter the adjustment.
%
%   ATTENTION:
%       -   The station and source coordinates are taken from the VieVS
%           internal structure arrays. They will match the TRF and CRF
%           coordinates (when these exist) rather than the coordinates in
%           the original NGS files.
%       -   The data quality flag is set to the value from the original NGS
%           file. It thus does not reflect the quality of the simulated
%           data (which would be always '0'). This means that if you
%           exclude observations with a quality flag higher than q,
%           the same observations will be eliminated in the original and
%           the simulated NGS files and the solution is run with exactly
%           the same number of observations.
%       -   Asterisks in the original NGS file will be substituted by
%           zeros in the simulated NGS file.
%
%   References:
%    ---
%
%   Input:
%      session   -    string with session name
%      dirpt     -    name of the subdirectory
%      yr        -    year
%      idays     -    number of days for which simulations were
%                     carried out = number of NGS files to be written
%      scan      -    the VieVS scan structure array
%      antenna   -    the VieVS antenna structure array
%      sources   -    the VieVS source structure array
%      sind      -    first index for the running number of the NGS file
%      zinp      -    zero input? (1 = yes, 0 = no) set in vie_sim.m
%      dirpt0    -    subdirectory in SIM
%   Output:
%      ---
%
%   External calls:
%      ---
%
%   Coded for VieVS:
%   June 2011 by Andrea Pany
%
%   Revision:
%   xx Oct 2015 by D. Mayer: added possibility to use subfolders
%   xx Nov 2016 by M. Schartner: new name of ngs file if it is a *_V* file
%   05 Dec 2018 by D. Landskron: clarification quality code / quality flag
%   09 May 2022 by L. Kern: simulations with vgosDB - bug fix
% ************************************************************************



function scan2ngs(session,antenna,scan,sources,yr,zinp,sind,nday,dirpt0)
tic
% check if SIM directory in DATA directory exists and create if not
if isempty(dirpt0)
    if ~isdir(['../DATA/SIM/', num2str(yr)])
        mkdir(['../DATA/SIM/' , num2str(yr)]);
    end
    sdir = ['../DATA/SIM/', num2str(yr)];
else
    if ~isdir(['../DATA/SIM/', dirpt0, '/',num2str(yr)])
        mkdir(['../DATA/SIM/', dirpt0, '/',num2str(yr)]);
    end
    sdir = ['../DATA/SIM/', dirpt0, '/',num2str(yr)];
end

for iday = 1:nday
    
    % output ngs file
    % for correct naming of the files (S001, S010, S100)
    % use the starting index specified by the user in the sim gui
    % if chara
    if length(session) > 9 % input = NGS file
        if strcmp(session(11),'V')
            if (iday+sind-1) < 10
                fn2 = strcat(session, '_S00', num2str(iday+sind-1));
            elseif (iday+sind-1) < 100
                fn2 = strcat(session, '_S0', num2str(iday+sind-1));
            else
                fn2 = strcat(session, '_S', num2str(iday+sind-1));
            end
        else
            if (iday+sind-1) < 10
                fn2 = strcat(session(1:10), 'S00', num2str(iday+sind-1));
            elseif (iday+sind-1) < 100
                fn2 = strcat(session(1:10), 'S0', num2str(iday+sind-1));
            else
                fn2 = strcat(session(1:10), 'S', num2str(iday+sind-1));
            end
        end
    else % input = vgosDB file
        if (iday+sind-1) < 10
            fn2 = strcat(session(1:9), '_S00', num2str(iday+sind-1));
        elseif (iday+sind-1) < 100
            fn2 = strcat(session(1:9), '_S0', num2str(iday+sind-1));
        else
            fn2 = strcat(session(1:9), '_S', num2str(iday+sind-1));
        end
    end
    fid2 = fopen([sdir,'/',fn2], 'w');
    
    % write header
    fprintf(fid2, '%s\n', ['DATA IN NGS FORMAT FROM DATABASE ',session]);
    fprintf(fid2, '%s\n', 'Simulated delays and observed rates in card #2');
    for i = 1:length(antenna)
        
        % station coordinates
        xx = sprintf('%15.5f',antenna(i).x);
        yy = sprintf('%15.5f',antenna(i).y);
        zz = sprintf('%15.5f',antenna(i).z);
        % axis offset
        oo = sprintf('%10.5f',antenna(i).offs);
        % write data to file
        if strcmp(antenna(i).axtyp,'X-Y1')
            axtyp = 'X-YN';
        else
            axtyp = antenna(i).axtyp;
        end
        line = ([antenna(i).name,'  ',xx,yy,zz,' ',axtyp,oo]);
        fprintf(fid2,'%s\n',line);
    end
    fprintf(fid2,'%s\n','$END');
    
    % write radio source position card
    for i = 1:length(sources)
        % convert to h/min/sec or deg/min/sec respectively
        [dra mra sra] = rad2houminsec(sources(i).ra2000);
        [dde mde sde] = rad2degminsec(sources(i).de2000);
        if (sources(i).de2000 < 0.0)   %%% 201212 12 Jing SUN%if (dde < 0.0)
            sign = '-';
        else
            sign = ' ';
        end
        adde = abs(dde);
        if (adde <= 9)
            sdde = ['0', num2str(adde)];
        else
            sdde = num2str(adde);
        end
        strde = [sign, sdde];
        % write data to file
        line = ([sources(i).name,'  ',sprintf('%2.0f',dra),' ',...
            sprintf('%2.0f',mra),' ',sprintf('%12.6f',sra),' ',...
            sprintf('%s',strde),' ',sprintf('%2.0f',mde),' ',...
            sprintf('%12.6f',sde)]);
        fprintf(fid2,'%s\n',line);
    end
    fprintf(fid2,'%s\n','$END');
    
    % write auxiliary parameters
    fprintf(fid2,'%s\n','$END');
    
    % write data cards
    seqnum = 0;
    for iscan = 1:length(scan)
        for iobs = 1:scan(iscan).nobs
            % sequence number
            seqnum = seqnum + 1;
            % observation epoch
            timyr = sprintf('%4.0f',scan(iscan).tim(1));
            if scan(iscan).tim(2) >= 10
                timmo = sprintf('%2.0f',scan(iscan).tim(2));
            else
                timmo = ['0',sprintf('%1.0f',scan(iscan).tim(2))];
            end
            if scan(iscan).tim(3) >= 10
                timda = sprintf('%2.0f',scan(iscan).tim(3));
            else
                timda = ['0',sprintf('%1.0f',scan(iscan).tim(3))];
            end
            if scan(iscan).tim(4) >= 10
                timhr = sprintf('%2.0f',scan(iscan).tim(4));
            else
                timhr = ['0',sprintf('%1.0f',scan(iscan).tim(4))];
            end
            if scan(iscan).tim(5) >= 10
                timmi = sprintf('%2.0f',scan(iscan).tim(5));
            else
                timmi = ['0',sprintf('%1.0f',scan(iscan).tim(5))];
            end
            timse = sprintf('%14.10f',scan(iscan).tim(6));
            
            %------------------------
            % write data card # 1
            %------------------------
            line = ([antenna(scan(iscan).obs(iobs).i1).name,'  ',...   % name of first station
                antenna(scan(iscan).obs(iobs).i2).name,'  ',...        % name of second station
                sources(scan(iscan).iso).name,' ',...                  % name of source
                timyr,' ',...                                          % year of obs
                timmo,' ',...                                          % month
                timda,' ',...                                          % day
                timhr,' ',...                                          % hour
                timmi,' ',...                                          % min
                timse,'          ',...                                 % sec
                sprintf('%8.0f',seqnum),'01']);                        % sequence number&card number
            fprintf(fid2,'%s\n',line);
            
            %-----------------------------------------------------------
            % write data card # 2 (delay, delay rate and formal errors,
            % quality flag)
            %-------------------------------------------
            % remove delay-corrections applied in read_ngs
            corcab = scan(iscan).stat(scan(iscan).obs(iobs).i2).cab - ...  % cable cal of second station
                scan(iscan).stat(scan(iscan).obs(iobs).i1).cab;       % cable cal of first station
            % zero input or simulated input?
            if zinp == 1
                % if zero input -> o = c -> o-c = 0
                vo = 1.0d9 * scan(iscan).obs(iobs).com + scan(iscan).obs(iobs).delion - corcab;
            else
                vo = 1.0d9 * (scan(iscan).obs(iobs).com + scan(iscan).obs(iobs).obs(1,iday)) + ...
                    scan(iscan).obs(iobs).delion - corcab;
            end
            % the sigma of the simulated delay observable is set to
            % the value of the simulated thermal noise (ionospheric formal
            % error has to be taken into account)
            sigdel = sqrt((1.0d9 * scan(iscan).obs(iobs).sig)^2 - ...
                scan(iscan).obs(iobs).sgdion^2);
            % write data
            line = ([sprintf('%20.8f',vo),sprintf('%10.5f',sigdel),...      % delay and formal delay error
                sprintf('%20.10f',0),sprintf('%10.5f',0),...           % delay rate and formal error (set to 0)
                sprintf('%2.0f',scan(iscan).obs(iobs).q_flag),' ',...  % data quality flag is set to the original value
                '  ',' ','  ','I',' ',sprintf('%8.0f',seqnum),'02']);
            fprintf(fid2,'%s\n',line);
            
            %--------------------------------------------
            % write data card # 3 (fringing information)
            %--------------------------------------------
            % all values are set to 0
            line = ([sprintf('%10.5f',0),sprintf('%10.5f',0),sprintf('%10.5f',0),...
                ' ',sprintf('%9.5f',0),sprintf('%20.15f',0),sprintf('%#10.0f',0),...
                sprintf('%8.0f',seqnum),'03']);
            fprintf(fid2,'%s\n',line);
            
            %------------------------------------------------------
            % write data card # 4 (system and antenna temperatures)
            %------------------------------------------------------
            % all values are set to 0
            line = ([sprintf('%10.2f',0),sprintf('%5.1f',0),sprintf('%10.2f',0),sprintf('%5.1f',0),...
                sprintf('%10.2f',0),sprintf('%5.1f',0),sprintf('%10.2f',0),sprintf('%5.1f',0),...
                '          ',sprintf('%8.0f',seqnum),'04']);
            fprintf(fid2,'%s\n',line);
            
            %-------------------------------------------
            % write data card # 5 (cable cal, WVR info)
            %-------------------------------------------
            % cable cal is taken from scan, all other parameters are set to 0
            line = ([sprintf('%10.5f',scan(iscan).stat(scan(iscan).obs(iobs).i1).cab), ... % cable cal station 1
                sprintf('%10.5f',scan(iscan).stat(scan(iscan).obs(iobs).i2).cab), ... % cable cal station 2
                sprintf('%10.5f',0),sprintf('%10.5f',0),sprintf('%10.5f',0),sprintf('%10.5f',0),...
                ' ',' ',' ',' ','      ',sprintf('%8.0f',seqnum),'05']);
            fprintf(fid2,'%s\n',line);
            
            %-------------------------------------------
            % write data card # 6 (met data)
            %-------------------------------------------
            line = ([sprintf('%10.3f',scan(iscan).stat(scan(iscan).obs(iobs).i1).temp),... % temperature station 1
                sprintf('%10.3f',scan(iscan).stat(scan(iscan).obs(iobs).i2).temp),... % temperature station 2
                sprintf('%10.3f',scan(iscan).stat(scan(iscan).obs(iobs).i1).pres),... % pressure station 1
                sprintf('%10.3f',scan(iscan).stat(scan(iscan).obs(iobs).i2).pres),... % pressure station 2
                sprintf('%10.3f',0),sprintf('%10.3f',0),' ','0',' ','0','      ',...
                sprintf('%8.0f',seqnum),'06']);
            fprintf(fid2,'%s\n',line);
            
            %-------------------------------------------
            % write data card # 8 (ionosphere corrections)
            %-------------------------------------------
            line = ([sprintf('%20.10f',scan(iscan).obs(iobs).delion),...
                sprintf('%10.5f',scan(iscan).obs(iobs).sgdion),...    % iono correction
                sprintf('%20.10f',0),sprintf('%10.5f',0),' ',...
                sprintf('%2.0f',scan(iscan).obs(iobs).q_flag_ion),'       ',...
                sprintf('%8.0f',seqnum),'08']);
            fprintf(fid2,'%s\n',line);
        end
    end
    fclose(fid2);
end

toc