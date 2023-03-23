% #########################################################################
% #     write_eopivs
% #########################################################################
%
% DESCRIPTION
%   Write EOP files in IVS format - *.eopi or *.eoxy
% 
%
% CREATED  
%   2021/05     Sigrid Boehm
%
% REFERENCES
%
%   External calls:
%   VieVS/COMMON/time/tai_utc.m
%   VieVS/COMMON/eop/rg_zont2.m
%	VieVS/VIE_SETUP/eop_get_approx.m
%
% INPUT
% - process_list            : Process list containing all sessions results of which should be written to the output file 
% - subdir                  : Sub-directory in VieVS file structure
% - outfile (optional)      : Output filename and directory
% - flag_intensive          : EOP from 24h or intensive sessions
% - flag_offsrate           : Output mode of daily estimated EOP (offset+rate or piecewise linear offsets)
%
% OUTPUT
%	creates txt-files in IVS EOP format with respect to input option
%   *.eopi file for intensive sessions, *.eoxy for 24h sessions
%
% CHANGES
% 04-2022 LK - implemented new IVS EOP format v3.0
% 11-2022 LK - add comment if estimated EOP are used in automated processing
% 12-2022 LK - implement changes of vgosDB naming convention + master file

function write_eopivs(process_list, subdir, outfile, flag_intensive, flag_offsrate, flag_incloutl)

    %% Get file name from GUI input
    % if a mat file is given instead of already loaded process list
    if strcmpi(process_list(end-3:end), '.mat')
	    temp=load(process_list);
	    fn=fieldnames(temp);
	    process_list=temp.(fn{1});
    end
    
    pl=size(process_list);
    
    if ~isempty(subdir)
	    if strcmpi(subdir(end), '/') || strcmpi(subdir(end), '\')
		    subdir=subdir(1:end-1);
	    end
    end
    
    if strcmp(outfile(end-3:end), '.txt')
        outfile = outfile(1:end-4);
    end
    
    %% Start collecting the data depending on session type (24h or intensive)
    
    % load superstationfile to get two character station code 
    load('../TRF/superstation.mat','superstations');
    code_name = [{superstations.code}' {superstations.name}'];
    clear superstations;
    curline = 0; 
    
    for j = 1:pl(1)
        
        % Get session name:
        % Unix:
        ind_unix = strfind(process_list(j,:), '/');
        % DOS:
        ind_dos = strfind(process_list(j,:), '\');
        ind = max([ind_unix, ind_dos]);
        sname = process_list(j,ind+1:end);
        sname = strtrim(sname); % Remove leading and trailing blanks from string - just to be sure...

        % Remove tags ([VGOSDB] or [VSO]) from string:
        char_id = strfind(sname, ' ');
        if ~isempty(char_id)
            sname = sname(1:char_id-1);
        end

        try 
            % load files
	        load(strcat('../DATA/LEVEL3/',subdir,'/x_',sname),'x_');
	        load(strcat('../DATA/LEVEL3/',subdir,'/opt_',sname),'opt_');
            load(strcat('../DATA/LEVEL3/',subdir,'/',sname,'_parameter'),'parameter');
        catch
            fprintf('Error: Check if x_, opt_, and _parameter.mat file exists for session %s!\n', sname)
            continue
        end

        % distinguish between old and new vgosDB naming convention
        if length(sname) == 15 && contains(sname, '-') % new vgosDB name YYYYMMDD-SSSSS
            syear = str2double(sname(1:4));
            scode = extractAfter(sname, '-');
        else % old vgosDB name YYMONDDSS
            % trim sname after files are loaded (remove NGS version number)
            if length(sname) > 9
                sname(10:end) = [];
            end
            syear = str2double(sname(1:2));
            if syear>70
                syear = num2str(1900+syear);
            else
                syear = num2str(2000+syear);
            end
            smondd = sname(3:7);
            sdbc = sname(8:9); sdbc = replace(sdbc,'_',' ');
        end
        
        %% get IVS session code from masterfile or vgosDB file
        path2masterfiles = '../DATA/MASTER/';
        
        if syear >= 2023
            IVSsesnam = scode;
        elseif flag_intensive && str2double(syear) > 1991
            path2masterfile = strcat(path2masterfiles, 'master', num2str(syear(3:4)), '-int.txt');

            % read intensive masterfile and save master struct
            master = readmasterfile(path2masterfile, syear, 'INT');

            % search for current session in masterfile using dbc and mondd
            flag_mondd = contains(master.mondd,smondd);
            flag_dbc = contains(master.dbc,sdbc);
            flag_sess = flag_mondd & flag_dbc;
            if sum(flag_sess) == 1 % session was found in the masterfile
                IVSsesnam = cell2mat(master.sessioncode(flag_sess));
            elseif sum(flag_sess) == 2 
                idx = find(flag_sess == 1);
                if contains(sname,'XU')
                    if strfind(master.sessioncode{idx(1)},'I')
                        flag_sess(idx(2)) = 0;
                    elseif strfind(master.sessioncode{idx(2)},'I')
                        flag_sess(idx(1)) = 0;
                    end
                    IVSsesnam = master.sessioncode{flag_sess};
                elseif contains(sname,'XT')
                    if strfind(master.sessioncode{idx(1)},'A')
                        flag_sess(idx(2)) = 0;
                    elseif strfind(master.sessioncode{idx(2)},'A')
                        flag_sess(idx(1)) = 0;
                    end
                    IVSsesnam = master.sessioncode{flag_sess};
                else
                    continue
                end
            elseif sum(flag_sess) == 0 % session was not found in intensive masterfile
                % search in masterfile for standard 24h sessions
                path2masterfile = strcat(path2masterfiles, 'master', syear(3:4), '.txt');
    
                % read standard masterfile and save master struct
                master = readmasterfile(path2masterfile, syear, '');
    
                % search for current session in masterfile using dbc and mondd
                flag_mondd = contains(master.mondd, smondd);
                if any(isspace(sdbc))
                    flag_dbc = strcmp(deblank(master.dbc),sdbc(1));
                else
                    flag_dbc = contains(master.dbc,sdbc);
                end
                flag_sess = flag_mondd & flag_dbc;
                if sum(flag_sess) == 1
                    IVSsesnam = master.sessioncode{flag_sess};
                else
                    % try VGOS masterfile - if one exists
                    try 
                        path2masterfile = strcat(path2masterfiles, 'master', syear(3:4), '-vgos.txt');
    
                        % read vgos masterfile and save master struct
                        master = readmasterfile(path2masterfile, syear, 'VGOS');
        
                        % search for current session in VGOS masterfile using dbc and mondd
                        flag_mondd = contains(master.mondd, smondd);
                        flag_dbc = contains(master.dbc, sdbc);
                        flag_sess = flag_mondd & flag_dbc;
                        if sum(flag_sess) == 1
                            IVSsesnam = master.sessioncode{flag_sess};
                        elseif sum(flag_sess) ~= 1 % session was not found in VGOS masterfile
                            IVSsesnam = checkvgosdb(sname, syear);
                        end
                    catch
                        IVSsesnam = checkvgosdb(sname, syear);
                    end
                end
            end
        else
            path2masterfile = strcat(path2masterfiles, 'master', syear(3:4), '.txt');

            % read standard masterfile and save master struct
            master = readmasterfile(path2masterfile, syear, '');
    
            % search for current session in masterfile using dbc and mondd
            flag_mondd = contains(master.mondd, smondd);
            if any(isspace(sdbc))
                flag_dbc = strcmp(deblank(master.dbc),sdbc(1));
            else
                flag_dbc = contains(master.dbc,sdbc);
            end
            flag_sess = flag_mondd & flag_dbc;
            if sum(flag_sess) == 1
                    IVSsesnam = master.sessioncode{flag_sess};
            elseif sum(flag_sess) ~= 1 % session was not found in masterfile
                % try VGOS masterfile - if one exists
                try 
                    path2masterfile = strcat(path2masterfiles, 'master', syear(3:4), '-vgos.txt');

                    % read vgos masterfile and save master struct
                    master = readmasterfile(path2masterfile, syear, 'VGOS');
    
                    % search for current session in VGOS masterfile using dbc and mondd
                    flag_mondd = contains(master.mondd, smondd);
                    flag_dbc = contains(master.dbc, sdbc);
                    flag_sess = flag_mondd & flag_dbc;
                    if sum(flag_sess) == 1
                        IVSsesnam = master.sessioncode{flag_sess};
                    elseif sum(flag_sess) ~= 1 % session was not found in VGOS masterfile
                        IVSsesnam = checkvgosdb(sname, syear);
                    end
                catch
                    IVSsesnam = checkvgosdb(sname, syear);
                end
            end
        end
        
        %% collect session specific information
        % WRMS
        wrms = x_.wrms*100/2.99792458; %[ps]
        if ((wrms > 999) || isnan(wrms))
            continue
        end
    
        % number of observations
        num_obs = x_.nobs;
    
        % station network
        netID = '';
        for ista = 1:length(x_.antenna)
            indcod = strcmp(code_name(:,2),x_.antenna(ista).name);
            if isempty(cell2mat(code_name(indcod,1)))
                continue
            end
            if ista == 1
                netID = strtrim(code_name(indcod));
            else
                netID = strcat(netID,'-',strtrim(code_name(indcod)));
            end
        end
        
        % get session duration and mjd of mid-session and of estimates
        sdur = (opt_.last_scan - opt_.first_scan)*24;
        mjdsesmid = (opt_.first_scan + opt_.last_scan)/2; 
        
        %% read and interpolate EOP 
        if parameter.lsmopt.xpol.model==1 % polar motion
           pol_est = true;
           if opt_.xpol.coef > 1.0e-4
               flag_poltight = false;
           else
               flag_poltight = true;
           end
           mjdp = x_.xpol.mjd;
           dmjdp = mjdp(2)-mjdp(1);
    
           % for 24h interval keep only the two values near session midpoint
           if (dmjdp>=1 && length(mjdp)>=3)
               if mjdp(2)==mjdsesmid
                   mjdp(2)=[];
               else
                   mjdp(abs(mjdp-mjdsesmid)>1)=[];
               end
           elseif (dmjdp<1 && flag_poltight)
               mjdp(2:end-1)=[];
               dmjdp = mjdp(2)-mjdp(1);
           end
        else
	       mjdp = [];
           dmjdp = 0;
           pol_est = false;
           flag_poltight = false;
        end
    
        if parameter.lsmopt.dut1.model==1 % dUT1
           ut1_est = true;
           if opt_.dut1.coef > 1.0e-4
               flag_ut1tight = false;
           else
               flag_ut1tight = true;
           end
           mjdu = x_.dut1.mjd;
           dmjdu = mjdu(2)-mjdu(1);
    
           % for 24h interval keep only the two values near session midpoint
           if (dmjdu >=1 && length(mjdu)>=3)
               if mjdu(2)==mjdsesmid
                   mjdu(2)=[];
               else
                   mjdu(abs(mjdu-mjdsesmid)>1)=[];
               end
           elseif (flag_intensive && length(mjdu)==3) % sometimes also intensives have three values
               mjdu(abs(mjdu-mjdsesmid)>dmjdu)=[];
           elseif (dmjdu<1 && flag_ut1tight && ~flag_intensive)
               mjdu(2:end-1)=[];
               dmjdu = mjdu(2)-mjdu(1);
           end
        else
	       mjdu = [];
           dmjdu = 0;
           ut1_est = false;
           flag_ut1tight = false;
        end
    
        % nutation if estimated
        if parameter.lsmopt.nutdx.model==1 % nutation
           nut_est = true;
           if opt_.nutdx.coef > 1.0e-4
               flag_nuttight = false;
           else
               flag_nuttight = true;
           end
           mjdn = x_.nutdx.mjd;
           dmjdn = mjdn(2)-mjdn(1);
    
           % for 24h interval keep only the two values near session midpoint
           if (dmjdn>=1 && length(mjdn)>=3)
               if mjdn(2)==mjdsesmid
                   mjdn(2)=[];
               else
                   mjdn(abs(mjdn-mjdsesmid)>1)=[];
               end
           elseif (dmjdn<1 && flag_nuttight)
               mjdn(2:end-1)=[];
               dmjdn = mjdn(2)-mjdn(1);
           end
        else
	       mjdn = [];
           dmjdn = 0;
           nut_est = false;
           flag_nuttight = false; 
        end
        
	    % make one vector with all mjd values
	    mjd = [mjdp mjdu mjdn];
	    mjd = sort(mjd);
	    mjd(diff(mjd)==0)=[];
        
	    % get common mjds
        indp = ismember(mjd , mjdp);
        indu = ismember(mjd , mjdu);
        indn = ismember(mjd , mjdn);
       
        if pol_est
           xp   = x_.xpol.val(ismember(x_.xpol.mjd,mjdp))';  
	       xp_e = x_.xpol.mx(ismember(x_.xpol.mjd,mjdp))';  
	       yp   = x_.ypol.val(ismember(x_.ypol.mjd,mjdp))';  
	       yp_e = x_.ypol.mx(ismember(x_.ypol.mjd,mjdp))';  
           if ~flag_incloutl
           % skip sessions if formal error/estimate is above threshold (10mas)
           if any(xp_e>10) || any(yp_e>10) || any(abs(xp)>10) || any(abs(yp)>10)
               continue
           end
           end
        end
        
        if ut1_est
	       ut   = x_.dut1.val(ismember(x_.dut1.mjd,mjdu))';   
	       ut_e = x_.dut1.mx(ismember(x_.dut1.mjd,mjdu))';  
           if ~flag_incloutl
           % skip sessions if formal error/estimate is above threshold (0.67ms)
           if any(ut_e>0.67) || any(abs(ut)>0.67) 
               continue
           end
           end
        end
        
        if nut_est
 	       dX   = x_.nutdx.val(ismember(x_.nutdx.mjd,mjdn))'; 
	       dX_e = x_.nutdx.mx(ismember(x_.nutdx.mjd,mjdn))'; 
	       dY   = x_.nutdy.val(ismember(x_.nutdy.mjd,mjdn))'; 
	       dY_e = x_.nutdy.mx(ismember(x_.nutdy.mjd,mjdn))'; 
           if ~flag_incloutl
           % skip sessions if formal error/estimate is above threshold (1.5mas)
           if any(dX_e>1.5) || any(dY_e>1.5) || any(abs(dX)>1.5) || any(abs(dY)>1.5)
               continue
           end
           end
        end 
        
	    % get a priori values for estimation epochs
	    [leapsec,leapepo]  = tai_utc(mjd,1);       % get difference TAI-UTC [s]
	    TT = mjd + (32.184 + leapsec')./86400;     % [MJD TT]
    
        [XP,YP,DUT1,DXap,DYap,~] = eop_get_approx(parameter,mjd,TT);
    
        % calculate total EOP values (tot = estimated + a priori)
        if ut1_est
            uttot = ut + DUT1(indu); % [ms]
            % Check if there was a leap second within the dates of the a priori
            % EOP. In load_eop the UT1-UTC values are homogenized for interpolation 
            % if there is a step.
            % This artificial leap seconds are corrected again before utfin is
            % written to the EOP file.
            MJDleapcheck =  parameter.eop.mjd;
            [~,leapepo]  = tai_utc(MJDleapcheck,1); 
            [~,indleap] = intersect(MJDleapcheck, leapepo);
            if ~isempty(intersect(MJDleapcheck, leapepo))
                isleap = true;
                switchleap = zeros(length(MJDleapcheck),1);
                if (indleap-1) < 6
                    switchleap(1:indleap-1)=-1;
                else
                    switchleap(indleap:end)=1;
                end
            else
                isleap = false;
            end
        end
        clear leapsec leapepo TT
    
        if pol_est
            xptot = xp + XP(indp);   % [mas]
            yptot = yp + YP(indp);   % [mas]
        end
            
        if nut_est
            dXtot = dX + DXap(indn); % [mas]
            dYtot = dY + DYap(indn); % [mas]
        end
        
        % for estimation intervals daily or longer and choice offset+rate the value is interpolated
        % linearly to the session middle and the corresponding rate is calculated
        if (dmjdp>=1 && flag_offsrate) || flag_poltight
            mjdp2_mid = mjdp(2)-mjdsesmid;
            mid_mjdp1 = mjdsesmid-mjdp(1);
            xpfin = xptot(1)*mjdp2_mid/dmjdp + xptot(2)*mid_mjdp1/dmjdp;
            ypfin = yptot(1)*mjdp2_mid/dmjdp + yptot(2)*mid_mjdp1/dmjdp; 
            if ~flag_incloutl
            if (any(abs(xpfin)>=0.65e3) || any(abs(ypfin)>=0.65e3))  % skip session if value is >0.65 as/s
                continue
            end
            end
            xpfin_e = sqrt((mjdp2_mid/dmjdp*xp_e(1))^2 + (mid_mjdp1/dmjdp*xp_e(2))^2); 
            ypfin_e = sqrt((mjdp2_mid/dmjdp*yp_e(1))^2 + (mid_mjdp1/dmjdp*yp_e(2))^2); 
            xprat = (xptot(2)-xptot(1))/dmjdp;
            yprat = (yptot(2)-yptot(1))/dmjdp;
            xprat_e = sqrt((xp_e(1)/dmjdp)^2 + (xp_e(2)/dmjdp)^2);
            yprat_e = sqrt((yp_e(1)/dmjdp)^2 + (yp_e(2)/dmjdp)^2);
            if flag_poltight
                xprat = NaN; yprat = NaN; xprat_e = NaN; yprat_e = NaN;
            end
            clear mjdp
            mjdp = mjdsesmid;
        elseif dmjdp == 0
            xpfin = NaN; ypfin = NaN; xpfin_e = NaN; ypfin_e = NaN;
            xprat = NaN; yprat = NaN; xprat_e = NaN; yprat_e = NaN;        
        elseif dmjdp>0 || ~flag_offsrate
            xpfin = xptot; ypfin = yptot; xpfin_e = xp_e; ypfin_e = yp_e;
            xprat = NaN; yprat = NaN; xprat_e = NaN; yprat_e = NaN;
            if ~flag_incloutl
            if (any(abs(xpfin)>=0.65e3) || any(abs(ypfin)>=0.65e3))  % skip session if value is >0.65 as/s
                continue
            end
            end
        end
        
        if (dmjdn>=1 && flag_offsrate) || flag_nuttight
            mjdn2_mid = mjdn(2)-mjdsesmid;
            mid_mjdn1 = mjdsesmid-mjdn(1);
            dXfin = dXtot(1)*mjdn2_mid/dmjdn + dXtot(2)*mid_mjdn1/dmjdn;
            dYfin = dYtot(1)*mjdn2_mid/dmjdn + dYtot(2)*mid_mjdn1/dmjdn; 
            if ~flag_incloutl
            if (any(abs(dXfin)>=5) || any(abs(dYfin)>=5))  % skip session if value is >5 mas
                continue
            end
            end
            dXfin_e = sqrt((mjdn2_mid/dmjdn*dX_e(1))^2 + (mid_mjdn1/dmjdn*dX_e(2))^2); 
            dYfin_e = sqrt((mjdn2_mid/dmjdn*dY_e(1))^2 + (mid_mjdn1/dmjdn*dY_e(2))^2); 
            dXrat = (dXtot(2)-dXtot(1))/dmjdn;
            dYrat = (dYtot(2)-dYtot(1))/dmjdn;
            dXrat_e = sqrt((dX_e(1)/dmjdn)^2 + (dX_e(2)/dmjdn)^2);
            dYrat_e = sqrt((dY_e(1)/dmjdn)^2 + (dY_e(2)/dmjdn)^2);
            if flag_nuttight
                dXrat = NaN; dYrat = NaN; dXrat_e = NaN; dYrat_e = NaN; 
            end
            clear mjdn
            mjdn = mjdsesmid;
        elseif dmjdn == 0
            dXfin = NaN; dYfin = NaN; dXfin_e = NaN; dYfin_e = NaN;
            dXrat = NaN; dYrat = NaN; dXrat_e = NaN; dYrat_e = NaN; 
        elseif dmjdn>0 || ~flag_offsrate
            dXfin = dXtot; dYfin = dYtot; dXfin_e = dX_e; dYfin_e = dY_e;
            dXrat = NaN; dYrat = NaN; dXrat_e = NaN; dYrat_e = NaN; 
            if ~flag_incloutl
            if (any(abs(dXfin)>=5) || any(abs(dYfin)>=5))  % skip session if value is >5 mas
                continue
            end
            end
        end
        
        if (dmjdu>=1 && flag_offsrate) || flag_intensive || flag_ut1tight
            mjdu2_mid = mjdu(2)-mjdsesmid;
            mid_mjdu1 = mjdsesmid-mjdu(1);
            [leapsec,~]  = tai_utc([mjdu(1) mjdu(2) mjdsesmid],1);       % get difference TAI-UTC [s]
	        TT = [mjdu(1) mjdu(2) mjdsesmid] + (32.184 + leapsec')./86400;     % [MJD TT]
            UT1corr  = rg_zont2(TT,2);  % [s]
            uttoti   = uttot - UT1corr(1:2,:)'*1e3;  % [ms]
            utfini = uttoti(1)*mjdu2_mid/dmjdu + uttoti(2)*mid_mjdu1/dmjdu;
            utfin = utfini + UT1corr(3,1)*1e3; % [ms]
            if ~flag_incloutl
            if any(abs(utfin)>=1e3) % skip session if value is >1 as/s
                continue
            end
            end
            utfin_e = sqrt((mjdu2_mid/dmjdu*ut_e(1))^2 + (mid_mjdu1/dmjdu*ut_e(2))^2); 
            utrat = (uttot(2)-uttot(1))/dmjdu;  
            lod = -utrat;
            utrat_e = sqrt((ut_e(1)/dmjdu)^2 + (ut_e(2)/dmjdu)^2);
            if flag_ut1tight
                lod = NaN; utrat = NaN; utrat_e = NaN; 
            end
            clear mjdu
            mjdu = mjdsesmid;
        elseif dmjdu == 0
            utfin = NaN; utfin_e = NaN; lod = NaN; utrat = NaN; utrat_e = NaN; 
        elseif dmjdu>0 || ~flag_offsrate
            utfin = uttot; utfin_e = ut_e; lod = NaN; utrat = NaN; utrat_e = NaN; 
            if ~flag_incloutl
            if any(abs(utfin)>=1e3) % skip session if value is >1 as/s
                continue
            end
            end
            if isleap
                [~,mjduind] = intersect(MJDleapcheck,mjdu);
                utfin = utfin + switchleap(mjduind)'*1e3;
            end
        end
            
        clear mjd indp indu indn
	    % make again one vector with all mjd values
	    mjd = [mjdp mjdu mjdn];
	    mjd = sort(mjd);
	    mjd(diff(mjd)==0)=[];
        
        % get common mjds
        indp = ismember(mjd , mjdp);
        indu = ismember(mjd , mjdu);
        indn = ismember(mjd , mjdn);
    
        filler = 'NA'; % parameter not estimated - place filler
        
        for mjdind = 1:length(mjd)
           eop{1}{curline+mjdind} = sprintf('%12.6f',mjd(mjdind)); % epoch [MJD]
           if indp(mjdind)
               if ~isnan(xpfin(1))
                    eop{2}{curline+mjdind}  = adjfmt(xpfin(mjdp==mjd(mjdind))*1e-3,'%9.7f'); % xPol [as]
                    eop{3}{curline+mjdind}  = adjfmt(ypfin(mjdp==mjd(mjdind))*1e-3,'%9.7f'); % yPol [as]
                    eop{7}{curline+mjdind}  = adjfmt(xpfin_e(mjdp==mjd(mjdind))*1e-3,'%9.7f'); % sig_xP [as]
                    eop{8}{curline+mjdind}  = adjfmt(ypfin_e(mjdp==mjd(mjdind))*1e-3,'%9.7f'); % sig_yP [as]
               else
                    eop{2}{curline+mjdind} = strjust(sprintf('%9s',filler),'center');
                    eop{3}{curline+mjdind} = strjust(sprintf('%9s',filler),'center');
                    eop{7}{curline+mjdind} = strjust(sprintf('%9s',filler),'center');
                    eop{8}{curline+mjdind} = strjust(sprintf('%9s',filler),'center');
               end               
               if ~isnan(xprat(1)) 
                    eop{20}{curline+mjdind} = adjfmt(xprat(mjdp==mjd(mjdind))*1e-3,'%10.8f'); % xPolR [as/d]
                    eop{21}{curline+mjdind} = adjfmt(yprat(mjdp==mjd(mjdind))*1e-3,'%10.8f'); % yPolR [as/d]
                    eop{25}{curline+mjdind} = adjfmt(xprat_e(mjdp==mjd(mjdind))*1e-3,'%10.8f'); % sig_xPR [as/d]
                    eop{26}{curline+mjdind} = adjfmt(yprat_e(mjdp==mjd(mjdind))*1e-3,'%10.8f'); % sig_yPR [as/d] 
               else
                    eop{20}{curline+mjdind} = strjust(sprintf('%10s',filler),'center');
                    eop{21}{curline+mjdind} = strjust(sprintf('%10s',filler),'center'); 
                    eop{25}{curline+mjdind} = strjust(sprintf('%10s',filler),'center');
                    eop{26}{curline+mjdind} = strjust(sprintf('%10s',filler),'center');            
                end
           else
                eop{2}{curline+mjdind} = strjust(sprintf('%9s',filler),'center');
                eop{3}{curline+mjdind} = strjust(sprintf('%9s',filler),'center');
                eop{7}{curline+mjdind} = strjust(sprintf('%9s',filler),'center');
                eop{8}{curline+mjdind} = strjust(sprintf('%9s',filler),'center');
                eop{20}{curline+mjdind} = strjust(sprintf('%10s',filler),'center');
                eop{21}{curline+mjdind} = strjust(sprintf('%10s',filler),'center'); 
                eop{25}{curline+mjdind} = strjust(sprintf('%10s',filler),'center');
                eop{26}{curline+mjdind} = strjust(sprintf('%10s',filler),'center');
           end
           if indu(mjdind)
               if ~isnan(utfin(1))
                    eop{4}{curline+mjdind} = adjfmt(utfin(mjdu==mjd(mjdind))*1e-3,'%10.8f'); % dUT1 [s]
                    eop{9}{curline+mjdind} = adjfmt(utfin_e(mjdu==mjd(mjdind))*1e-3,'%10.8f'); % sig_UT [s]
               else
                    eop{4}{curline+mjdind} = strjust(sprintf('%10s',filler),'center');
                    eop{9}{curline+mjdind} = strjust(sprintf('%10s',filler),'center');
               end
               if ~isnan(lod(1))  
                    eop{22}{curline+mjdind} = adjfmt(lod(mjdu==mjd(mjdind))*1e-3,'%11.9f'); % LOD [s/d]
                    eop{27}{curline+mjdind} = adjfmt(utrat_e(mjdu==mjd(mjdind))*1e-3,'%11.9f'); % sig_LOD [s/d] 
               else
                    eop{22}{curline+mjdind} = strjust(sprintf('%11s',filler),'center');
                    eop{27}{curline+mjdind} = strjust(sprintf('%11s',filler),'center');     
               end
           else
               eop{4}{curline+mjdind} = strjust(sprintf('%10s',filler),'center');
               eop{9}{curline+mjdind} = strjust(sprintf('%10s',filler),'center');
               eop{22}{curline+mjdind} = strjust(sprintf('%11s',filler),'center');
               eop{27}{curline+mjdind} = strjust(sprintf('%11s',filler),'center');        
           end
           if indn(mjdind)
               if ~isnan(dXfin(1))
                    eop{5}{curline+mjdind} = adjfmt(dXfin(mjdn==mjd(mjdind)),'%6.4f'); % dPsi [mas]
                    eop{6}{curline+mjdind} = adjfmt(dYfin(mjdn==mjd(mjdind)),'%6.4f'); % dEps [mas]
                    eop{10}{curline+mjdind} = adjfmt(dXfin_e(mjdn==mjd(mjdind)),'%8.4f'); % sig_dPsi [mas]
                    eop{11}{curline+mjdind} = adjfmt(dYfin_e(mjdn==mjd(mjdind)),'%8.4f'); % sig_dEps [mas]
               else 
                    eop{5}{curline+mjdind} = strjust(sprintf('%6s',filler),'center');
                    eop{6}{curline+mjdind} = strjust(sprintf('%6s',filler),'center');
                    eop{10}{curline+mjdind} = strjust(sprintf('%8s',filler),'center');
                    eop{11}{curline+mjdind} = strjust(sprintf('%8s',filler),'center');
               end
               if ~isnan(dXrat(1))
                    eop{23}{curline+mjdind} = adjfmt(dXrat(mjdn==mjd(mjdind)),'%7.5f'); % dPsiR [mas/d]
                    eop{24}{curline+mjdind} = adjfmt(dYrat(mjdn==mjd(mjdind)),'%7.5f'); % dEpsR [mas/d]
                    eop{28}{curline+mjdind} = adjfmt(dXrat_e(mjdn==mjd(mjdind)),'%7.5f'); % sig_dPsiR [mas/d]
                    eop{29}{curline+mjdind} = adjfmt(dYrat_e(mjdn==mjd(mjdind)),'%7.5f'); % sig_dEpsR [mas/d] 
               else
                   eop{23}{curline+mjdind} = strjust(sprintf('%7s',filler),'center');
                   eop{24}{curline+mjdind} = strjust(sprintf('%7s',filler),'center');
                   eop{28}{curline+mjdind} = strjust(sprintf('%7s',filler),'center');
                   eop{29}{curline+mjdind} = strjust(sprintf('%7s',filler),'center');
               end
           else
               eop{5}{curline+mjdind} = strjust(sprintf('%6s',filler),'center');
               eop{6}{curline+mjdind} = strjust(sprintf('%6s',filler),'center');
               eop{10}{curline+mjdind} = strjust(sprintf('%8s',filler),'center');
               eop{11}{curline+mjdind} = strjust(sprintf('%8s',filler),'center');
               eop{23}{curline+mjdind} = strjust(sprintf('%7s',filler),'center');
               eop{24}{curline+mjdind} = strjust(sprintf('%7s',filler),'center');
               eop{28}{curline+mjdind} = strjust(sprintf('%7s',filler),'center');
               eop{29}{curline+mjdind} = strjust(sprintf('%7s',filler),'center');
           end
           eop{12}{curline+mjdind} = sprintf('%8.2f',wrms); % WRMS [ps]
           eop{13}{curline+mjdind} = strjust(sprintf('%9s',filler),'center'); % cor_xPyP [-]
           eop{14}{curline+mjdind} = strjust(sprintf('%9s',filler),'center'); % cor_xPUT [-]
           eop{15}{curline+mjdind} = strjust(sprintf('%9s',filler),'center'); % cor_xPUT [-]
           eop{16}{curline+mjdind} = strjust(sprintf('%9s',filler),'center'); % cor_dPdE [-]
           eop{17}{curline+mjdind} = sprintf('%6.0f',num_obs); % Nobs [-]
           eop{18}{curline+mjdind} = sprintf('%8s',IVSsesnam); % SessID [-]
           eop{19}{curline+mjdind} = sprintf('%7.2f',sdur); % Span [h]
           eop{30}{curline+mjdind} = sprintf('%-64s',string(netID)); % Network [-]
           eop{31}{curline+mjdind} = sprintf('%-64s',''); % Comment [-]
           
        end        
        curline = curline+mjdind;
    end
    
    % sort output cell array
    [~,idx] = sort(eop{1});
    for n = 1:size(eop,2)
         eop{n} = eop{n}(idx);
    end
    
    if flag_intensive
         eopout_name=[outfile,'.eopi'];
    else
         eopout_name=[outfile,'.eoxy'];
    end
    
    %% generating EOP file
    fprintf('Write IVS EOP file ...\n');
    generation_time = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')); 
    data_start = char(datetime(str2double(eop{1}{1}), 'ConvertFrom', 'modifiedjuliandate', 'Format','yyyy-MM-dd HH:mm:ss'));
    data_end = char(datetime(str2double(eop{1}{end}), 'ConvertFrom', 'modifiedjuliandate', 'Format','yyyy-MM-dd HH:mm:ss'));
    num_entries = size(eop{1},2);
    fid = fopen(eopout_name, 'w');
    
    % data description line
    fprintf(fid,'%s %sT%s %s %sT%s %sT%s %s\n','%=IVS-EOP 3.0 VIE',generation_time(1:10),generation_time(12:end),'VIE',data_start(1:10),data_start(12:end),data_end(1:10),data_end(12:end),'UTC R');
    
    % header
    fprintf(fid,'%s\n','+HEADER');
    fprintf(fid,'%s\t\t%sT%s\n','GENERATION_TIME', generation_time(1:10),generation_time(12:end));
    fprintf(fid,'%s\t\t%sT%s\n','DATA_START', data_start(1:10),data_start(12:end));
    fprintf(fid,'%s\t\t%sT%s\n','DATA_END', data_end(1:10),data_end(12:end));
    fprintf(fid,'%s\t\t%s\n','DESCRIPTION', 'Earth orientation parameters from the VLBI single session solution');
    fprintf(fid,'%s\t\t%s\n','ANALYSIS_CENTER', 'NA');
    fprintf(fid,'%s\t\t\t%s\n','CONTACT', 'NA');
    fprintf(fid,'%s\t\t%s\n','SOFTWARE', 'VieVS v3.2');
    
    % session type - technique
    if flag_intensive
        fprintf(fid,'%s\t\t%s\n','TECHNIQUE', 'VINT'); % VLBI INT 
    else
        fprintf(fid,'%s\t\t%s\n','TECHNIQUE', 'VLBI'); % VLBI R1/R4 + VGOS + others
    end
    fprintf(fid,'%s\t\t%s\n','NUTATION_TYPE', 'CIO-BASED ');
    fprintf(fid,'%s\t\t%s\n','ROTATION_TYPE', 'UT1-UTC_LOD');
    
    % reference frames - trf/crf_apriori
    fprintf(fid,'%s\t\t%s\n','TRF_APRIORI', parameter.vie_init.trf{2});
    fprintf(fid,'%s\t\t%s\n','CRF_APRIORI', parameter.vie_init.crf{2});
    
    % EOP apriori
    if strcmp(parameter.vie_mod.eopoc,'IERS_Desai_Sibois.dat')
        fprintf(fid,'%s\t\t%s\n','EOP_SUBDAILY', 'DESAI-SIBOIS');
    else
        fprintf(fid,'%s\t\t%s\n','EOP_SUBDAILY', parameter.vie_mod.eopoc);
    end
    if strcmp(parameter.vie_mod.EOPfile,'finals_all_IAU2000.txt')
        fprintf(fid,'%s\t\t%s\n','EOP_APRIORI', 'finals2000A.all');
    else
        fprintf(fid,'%s\t\t%s\n','EOP_APRIORI', parameter.vie_mod.EOPfile);
    end
    
    % EOP estimated
    param_deg = '_BSP_1'; % (BSP=BSpline, 1=Degree)
    if flag_intensive 
        if flag_ut1tight == true
            param_deg = ''; % offset
            fprintf(fid,'%s\t\t%s\t\t%5.3f %s\n','EOP_ESTIMATED', ['DUT1' param_deg], parameter.lsmopt.dut1.coef*1000, 'uas');
        else
            param_deg = '_BSP_1'; % (BSP=BSpline, 1=Degree)
            fprintf(fid,'%s\t\t%s\t%5.3f %s\n','EOP_ESTIMATED', ['DUT1' param_deg], parameter.lsmopt.dut1.coef*1000, 'uas');
        end
    else
        if parameter.lsmopt.xpol.model == 1
            if flag_poltight == true
                param_deg = ''; % offset
                fprintf(fid,'%s\t\t%s\t\t%5.3f %s\n','EOP_ESTIMATED', ['XPOL' param_deg], parameter.lsmopt.xpol.coef, 'mas');
            else
                param_deg = '_BSP_1'; % (BSP=BSpline, 1=Degree)
                fprintf(fid,'%s\t\t%s\t%5.3f %s\n','EOP_ESTIMATED', ['XPOL' param_deg], parameter.lsmopt.xpol.coef, 'mas');
            end
        else
            fprintf(fid,'%s\t\t%s\t\t%s %s\n','EOP_ESTIMATED', 'XPOL', 'NONE', 'mas');
        end
        if parameter.lsmopt.ypol.model == 1
            if flag_poltight == true
                param_deg = ''; % offset
                fprintf(fid,'%s\t\t%s\t\t%5.3f %s\n','EOP_ESTIMATED', ['YPOL' param_deg], parameter.lsmopt.ypol.coef, 'mas');
            else
                param_deg = '_BSP_1'; % (BSP=BSpline, 1=Degree)
                fprintf(fid,'%s\t\t%s\t%5.3f %s\n','EOP_ESTIMATED', ['YPOL' param_deg], parameter.lsmopt.ypol.coef, 'mas');
            end
        else
            fprintf(fid,'%s\t\t%s\t\t%s %s\n','EOP_ESTIMATED', 'YPOL', 'NONE', 'mas');
        end
        if parameter.lsmopt.dut1.model == 1
            if flag_ut1tight == true
                param_deg = ''; % offset
                fprintf(fid,'%s\t\t%s\t\t%5.3f %s\n','EOP_ESTIMATED', ['DUT1' param_deg], parameter.lsmopt.dut1.coef, 'mas');
            else
                param_deg = '_BSP_1'; % (BSP=BSpline, 1=Degree)
                fprintf(fid,'%s\t\t%s\t%5.3f %s\n','EOP_ESTIMATED', ['DUT1' param_deg], parameter.lsmopt.dut1.coef, 'mas');
            end
        else
            fprintf(fid,'%s\t\t%s\t\t%s %s\n','EOP_ESTIMATED', 'DUT1', 'NONE', 'mas');
        end
        if parameter.lsmopt.nutdx.model == 1
            if flag_nuttight == true
                param_deg = ''; % offset
                fprintf(fid,'%s\t\t%s\t\t%5.3f %s\n','EOP_ESTIMATED', ['DX' param_deg], parameter.lsmopt.nutdx.coef*1000, 'uas');
            else
                param_deg = '_BSP_1'; % (BSP=BSpline, 1=Degree)
                fprintf(fid,'%s\t\t%s\t%5.3f %s\n','EOP_ESTIMATED', ['DX' param_deg], parameter.lsmopt.nutdx.coef*1000, 'uas');
            end
        else
            fprintf(fid,'%s\t\t%s\t\t%s %s\n','EOP_ESTIMATED', ['DX' param_deg], 'NONE', 'uas');
        end
        if parameter.lsmopt.nutdy.model == 1
            if flag_nuttight == true
                param_deg = ''; % offset
                fprintf(fid,'%s\t\t%s\t\t%5.3f %s\n','EOP_ESTIMATED', ['DY' param_deg], parameter.lsmopt.nutdy.coef*1000, 'uas');
            else
                param_deg = '_BSP_1'; % (BSP=BSpline, 1=Degree)
                fprintf(fid,'%s\t\t%s\t%5.3f %s\n','EOP_ESTIMATED', ['DY' param_deg], parameter.lsmopt.nutdy.coef*1000, 'uas');
            end
        else
            fprintf(fid,'%s\t\t%s\t\t%s %s\n','EOP_ESTIMATED', ['DY' param_deg], 'NONE', 'uas');
        end
    end
    fprintf(fid,'%s\t%.0f\n','NUMBER_OF_ENTRIES', num_entries);			
    fprintf(fid,'%s\n','-HEADER');
    
    % data
    fprintf(fid,'%s\n','+DATA');
    fprintf(fid,'%s\n','# All fields are in free format separated by blanks.');          
    fprintf(fid,'%s\n','# If a parameter was not estimated a filler NA is placed.');
    fprintf(fid,'%s\n','#     1           2        3          4        5      6       7         8         9          10       11       12       13        14        15        16       17     18       19       20         21          22       23      24        25        26          27        28      29          30                                                                31');
    fprintf(fid,'%s\n','#   epoch       xPol      yPol      dUT1      dX     dY     sig_xP    sig_yP    sig_UT     sig_dX   sig_dY    wRMS   cor_xPyP  cor_xPUT  cor_yPUT  cor_dXdY   nObs   sessID   span      xPolR      yPolR       LOD      dXR     dYR    sig_xPR    sig_yPR     sig_LOD   sig_dXR sig_dYR     network                                                          comments');
    fprintf(fid,'%s\n','#   [MJD]       [as]      [as]       [s]     [mas]  [mas]    [as]      [as]      [s]        [mas]    [mas]    [ps]     [-]       [-]       [-]       [-]      [-]     [-]      [h]     [as/d]     [as/d]       [s]    [mas/d] [mas/d]   [as/d]     [as/d]       [s]     [mas/d] [mas/d]       [-]                                                               [-]');
    for n = 1:size(eop{1},2)
        fprintf(fid,'%12s %+9s %+9s %+10s %+6s %+6s %+9s %+9s %+10s %+8s %+8s %+8s %+9s %+9s %+9s %+9s %+6s %8s %+7s %+10s %+10s %+11s %+7s %+7s %+10s %+10s %+11s %+7s %+7s     %-64s %-64s\n', ...
                eop{1}{n},eop{2}{n},eop{3}{n},eop{4}{n},eop{5}{n},...
                eop{6}{n},eop{7}{n},eop{8}{n},eop{9}{n},eop{10}{n},...
                eop{11}{n},eop{12}{n},eop{13}{n},eop{14}{n},eop{15}{n},...
                eop{16}{n},eop{17}{n},eop{18}{n},eop{19}{n},eop{20}{n},...
                eop{21}{n},eop{22}{n},eop{23}{n},eop{24}{n},eop{25}{n},...
                eop{26}{n},eop{27}{n},eop{28}{n},eop{29}{n},eop{30}{n},eop{31}{n});
    end
    fprintf(fid,'%s\n', '-DATA');
    fprintf(fid,'%s\n', '%IVS-EOP 3.0 END');
    fclose(fid);
    fprintf('Finished!\n');
end

function [master] = readmasterfile(path2masterfile, year, type)
    % DESCRIPTION
    %   reads masterfiles
    %
    % INPUT
    % - path2masterfile  : path to masterfile
    % - year of session and type ('INT', 'VGOS', else: 24h SX session)
    %
    % OUTPUT
    % - master           : struct with sessioncode, mondd, dbc

    % read masterfile
    fid = fopen(path2masterfile);   
    master = struct();  

    if strcmp(type,'INT')
        % for sessions from < 2019 number of headerlines is 10 (for >= 2019 use 8)
        if str2double(year) < 2019
            data = textscan(fid,'|%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter', '|','Headerlines', 10, 'CommentStyle', '---');
        elseif str2double(year) >= 2019
            data = textscan(fid,'|%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter', '|','Headerlines', 8, 'CommentStyle', '---');
        end
    elseif strcmp(type,'VGOS')
        data = textscan(fid,'|%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter', '|','Headerlines', 8, 'CommentStyle', '---');
    else
        if str2double(year) > 2017 % masterfiles since 2018 have less headerlines
            data = textscan(fid,'|%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter', '|','Headerlines', 8, 'CommentStyle', '---');
        else 
            data = textscan(fid,'|%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter', '|','Headerlines', 9, 'CommentStyle', '---');
        end
    end
    
    master.sessioncode = data{2};
    master.mondd = data{3};
    master.dbc = data{12};
    fclose(fid);
end

function [IVSsesnam] = checkvgosdb(sname, syear)
% if session was not found in master file try to get IVS session code from vgosDB file
    if contains(sname,'_')
        ind_ul = strfind(sname,'_');
        sname(ind_ul:end)=[];
    end
    vgosdbfile = strcat('../DATA/vgosDB/',syear,'/',sname,'.tgz');
    if exist(vgosdbfile,'file')
        untar(vgosdbfile);
        try
            IVSsesnam = (ncread(strcat(sname,'/Head.nc'),'ExpName')');
        catch
            try
                IVSsesnam = (ncread(strcat(sname,'/Head.nc'),'Session')');
            catch
                fprintf('ERROR: IVS code for session %s is not available, UNAVBL is placed instead of IVS session code!\n',sname);
                IVSsesnam = 'UNAVBL';    
            end
        end
        rmdir(sname, 's');
    else
        fprintf('ERROR: IVS code for session %s is not available, UNAVBL is placed instead of IVS session code!\n',sname);
        IVSsesnam = 'UNAVBL';
    end
end

function [outstr] = adjfmt(value,fmt)
    % DESCRIPTION
    %   adjusts format and removes leading zero if negative (-0.5 -> -.5)
    %
    % INPUT
    % - value           : input value
    % - fmt             : output format (e.g. '3.5f')
    %
    % OUTPUT
    % - outsr           : output sring
    
    % check if value is negative
    if value < 0
        dig = extractAfter(extractBefore(fmt,'.'),'%'); % number of digits
        dec = str2num(extractBefore(extractAfter(fmt,'.'),'f')); % number of decimals
            
        if value <= -1 % reduce number of decimals
            dec = string(dec-1); % reduce number by one to make room for '-'
            fmt = strcat('%',dig,'.',dec,'f');
            outstr = sprintf(fmt,value);
    
        else % remove leading zero to make room for '-'
            outstr = sprintf(fmt,value);
            outstr(2) = []; 
        end
    else
        outstr = sprintf(fmt,value);
    end
end
