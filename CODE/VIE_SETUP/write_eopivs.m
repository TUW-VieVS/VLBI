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
%		/VIE_SETUP/eop_get_approx.m
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
function write_eopivs(process_list, subdir, outfile,flag_intensive,flag_offsrate)

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
filler = '-0';   

for j = 1:pl(1)
    
    % Get session name:
    % Unix:
    ind_unix = strfind(process_list(j,:), '/');
    % DOS:
    ind_dos = strfind(process_list(j,:), '\');
    ind = max([ind_unix, ind_dos]);
    sname = process_list(j,ind+1:end);
    sname = strtrim(sname); % Remove leading and trailing blanks from string - just to be sure...
    yy = str2double(sname(1:2));
    if yy>70
        syear = num2str(1900+yy);
    else
        syear = num2str(2000+yy);
    end
    
    % Remove tags ([VGOSDB] or [VSO]) from string:
    char_id = strfind(sname, ' ');
    if ~isempty(char_id)
        sname = sname(1:char_id-1);
    end
    
    smondd = sname(3:7);
    sdbc = sname(8:9);
    
	load(strcat('../DATA/LEVEL3/',subdir,'/x_',sname),'x_');
	load(strcat('../DATA/LEVEL3/',subdir,'/opt_',sname),'opt_');
    load(strcat('../DATA/LEVEL3/',subdir,'/',sname,'_parameter'),'parameter');
    
        %% get IVS session code from masterfile or vgosDB file
        filepath_masterfiles = '../DATA/MASTER/';
        
        if flag_intensive
            filepath_masterfile = strcat(filepath_masterfiles, 'master', num2str(yy,'%02d'), '-int.txt');

            % read masterfile
            fid = fopen(filepath_masterfile);   
            master = struct();  

            % for sessions from < 2019 number of headerlines is 10 (for >= 2019 use 8)
            if yy < 19
                data = textscan(fid,'|%s %s %s %s %f:%f %f %s %s %s %s |%s %s %f %s %s','Delimiter', '|','Headerlines', 10, 'CommentStyle', '---');
            elseif yy >= 19
                data = textscan(fid,'|%s %s %s %s %f:%f %f %s %s %s %s |%s %s %f %s %s','Delimiter', '|','Headerlines', 8, 'CommentStyle', '---');
            end

            master.sessioncode = data{2};
            master.mondd = data{3};
            master.dbc = data{12};
            fclose(fid);

            % search for selected session in masterfile using dbc and mondd
            flag_mondd = contains(master.mondd,smondd);
            flag_dbc = contains(master.dbc,sdbc);
            current_sess = flag_mondd & flag_dbc;
            if sum(current_sess) == 1 % session was found in the masterfile
                IVSsesnam = cell2mat(master.sessioncode(current_sess));
            end
            
        else
        
        % find right masterfile for selected session
        stryy = sprintf('%02d',yy);
        filepath_masterfile = strcat(filepath_masterfiles,'master', stryy,'.txt');

        % read masterfile 
        fid = fopen(filepath_masterfile);
        master = struct(); 
        
        if syear > 2017 % masterfiles since 2018 have less headerlines
            header = textscan(fid,'|%s %s %s %s %f:%f %f %s %s %s %s %f %s %s %f %s','Delimiter', '|','Headerlines', 8, 'CommentStyle', '---');
        else
            header = textscan(fid,'|%s %s %s %s %f:%f %f %s %s %s %s %f %s %s %f %s','Delimiter', '|','Headerlines', 9, 'CommentStyle', '---');
        end

        master.sessioncode = header{2};
        master.mondd = header{3};
        master.dbc = header{13};
        fclose(fid);

        % read VGOS masterfile (for some years seperate VGOS masterfiles exist)
        filepath_masterfile_vgos = strcat(filepath_masterfiles,'master', stryy,'-vgos.txt');
        vgosmaster = false;
        if exist(filepath_masterfile_vgos,'file') 
            vgosmaster = true;
            fid = fopen(filepath_masterfile_vgos);
            header = textscan(fid,'|%s %s %s %s %f:%f %f %s %s %s %s %f %s %s %f %s','Delimiter', '|','Headerlines', 8, 'CommentStyle', '---');
            masterVGOS = struct();  
            masterVGOS.sessioncode = header{2};
            masterVGOS.mondd = header{3};
            masterVGOS.dbc = header{13};
            fclose(fid);
        end

        % search for selected session in masterfile using dbc and mondd
        flag_mondd = contains(master.mondd,smondd);
        flag_dbc = contains(master.dbc,sdbc);
        current_sess = flag_mondd & flag_dbc;
        if sum(current_sess) == 1 % session was found in the masterfile
            IVSsesnam = cell2mat(master.sessioncode(current_sess));
        elseif vgosmaster % session was not found in the masterfile, check the VGOS masterfile 
            flag_mondd = contains(masterVGOS.mondd,smondd);
            flag_dbc = contains(masterVGOS.dbc,sdbc);
            current_sess = flag_mondd & flag_dbc;
            if sum(current_sess) == 1 
                IVSsesnam = cell2mat(masterVGOS.sessioncode(current_sess));
            end
        end
        end
        
        % if session was not found in master file try to get IVS session code from vgosDB file
        if sum(current_sess) ~= 1
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
    
    %% collect session specific information
    wrms_pfr = x_.wrms*100/2.99792458; %[ps]
    if wrms_pfr > 999
        continue
    end
    num_obs = x_.nobs;
    netID = '';
    for ista = 1:length(x_.antenna)
     indcod = strcmp(code_name(:,2),x_.antenna(ista).name);
     if isempty(cell2mat(code_name(indcod,1)))
         continue
     end
     netID = strcat(netID,strtrim(code_name(indcod)));
    end
    
    % get session duration and mjd of mid-session and of estimates
    ses_dur = (opt_.last_scan - opt_.first_scan)*24;
    mjdsesmid = (opt_.first_scan + opt_.last_scan)/2; 
    
    %% read and interpolate EOP 
    if parameter.lsmopt.xpol.model==1
       pol_est = true;
       if opt_.xpol.coef > 1.0e-4
           flag_poltight = false;
       else
           flag_poltight = true;
       end
       mjdp = x_.xpol.mjd;
       dmjdp = mjdp(2)-mjdp(1);
       % for 24h interval keep only the two values near session midpoint
       % if offset+rate is selected or parameter was estimated with tight constraints
       if (dmjdp>=1 && length(mjdp)==3 && (flag_offsrate || flag_poltight))
           mjdp(abs(mjdp-mjdsesmid)>1)=[];
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

	if parameter.lsmopt.dut1.model==1
       ut1_est = true;
       if opt_.dut1.coef > 1.0e-4
           flag_ut1tight = false;
       else
           flag_ut1tight = true;
       end
       mjdu = x_.dut1.mjd;
       dmjdu = mjdu(2)-mjdu(1);
       % for 24h interval keep only the two values near session midpoint
       % if offset+rate is selected or parameter was estimated with tight constraints
       if (dmjdu >=1 && length(mjdu)==3 && (flag_offsrate || flag_ut1tight))
           mjdu(abs(mjdu-mjdsesmid)>1)=[];
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
    if parameter.lsmopt.nutdx.model==1
       nut_est = true;
       if opt_.nutdx.coef > 1.0e-4
           flag_nuttight = false;
       else
           flag_nuttight = true;
       end
       mjdn   = x_.nutdx.mjd;
       dmjdn = mjdn(2)-mjdn(1);
       % for 24h interval keep only the two values near session midpoint
       % if offset+rate is selected or parameter was estimated with tight constraints
       if (dmjdn>=1 && length(mjdn)==3 && (flag_offsrate || flag_nuttight))
           mjdn(abs(mjdn-mjdsesmid)>1)=[];
       elseif (dmjdn<1 && flag_nuttight)
           mjdn(2:end-1)=[];
           dmjdn = mjdn(2)-mjdn(1);
       end
    else
	   mjdn   = [];
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
    end
    
	if ut1_est
	   ut   = x_.dut1.val(ismember(x_.dut1.mjd,mjdu))';   
	   ut_e = x_.dut1.mx(ismember(x_.dut1.mjd,mjdu))';  
    end
    
    if nut_est
 	   dX   = x_.nutdx.val(ismember(x_.nutdx.mjd,mjdn))'; 
	   dX_e = x_.nutdx.mx(ismember(x_.nutdx.mjd,mjdn))'; 
	   dY   = x_.nutdy.val(ismember(x_.nutdy.mjd,mjdn))'; 
	   dY_e = x_.nutdy.mx(ismember(x_.nutdy.mjd,mjdn))'; 
    end 
    
	% get a priori values for estimation epochs
	[leapsec,leapepo]  = tai_utc(mjd,1);       % get difference TAI-UTC [s]
	TT = mjd + (32.184 + leapsec')./86400;     % [MJD TT]

    [XP,YP,DUT1,DXap,DYap,~] = eop_get_approx(parameter,mjd,TT);
    
    % calculate total EOP values (tot = estimated + a priori)
    if ut1_est
    uttot = ut + DUT1(indu); % [ms]
    if ~isempty(intersect(leapepo,mjd)) && any(diff(leapsec)~=0)
        [~,~,IB]=intersect(leapepo,mjd);
    	uttot(IB:end) = ut(IB:end) + DUT1(IB:end) +1e3; % [ms]
    end
    end
    clear leapsec leapepo TT

    if pol_est
    xptot = xp + XP(indp);   % [mas]
    yptot = yp + YP(indp);   % [mas]
    end
        %
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
    end
    
    if (dmjdn>=1 && flag_offsrate) || flag_nuttight
        mjdn2_mid = mjdn(2)-mjdsesmid;
        mid_mjdn1 = mjdsesmid-mjdn(1);
        dXfin = dXtot(1)*mjdn2_mid/dmjdn + dXtot(2)*mid_mjdn1/dmjdn;
        dYfin = dYtot(1)*mjdn2_mid/dmjdn + dYtot(2)*mid_mjdn1/dmjdn; 
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
    
    for mjdind = 1:length(mjd)
       eopout_line{1}{curline+mjdind} = sprintf('%12.6f',mjd(mjdind));
       if indp(mjdind)
           if ~isnan(xpfin(1))
               eopout_line{2}{curline+mjdind} = sprintf('%9.6f',xpfin(mjdp==mjd(mjdind))*1e-3); % [as]
               eopout_line{3}{curline+mjdind} = sprintf('%9.6f',ypfin(mjdp==mjd(mjdind))*1e-3); % [as]
               eopout_line{7}{curline+mjdind} = sprintf('%9.6f',xpfin_e(mjdp==mjd(mjdind))*1e-3); % [as]
               eopout_line{8}{curline+mjdind} = sprintf('%9.6f',ypfin_e(mjdp==mjd(mjdind))*1e-3); % [as]
           else
               eopout_line{2}{curline+mjdind} = sprintf('%9s',filler);
               eopout_line{3}{curline+mjdind} = sprintf('%9s',filler);
               eopout_line{7}{curline+mjdind} = sprintf('%9s',filler);
               eopout_line{8}{curline+mjdind} = sprintf('%9s',filler);
           end               
           if ~isnan(xprat(1))
             eopout_line{20}{curline+mjdind} = sprintf('%9.6f',xprat(mjdp==mjd(mjdind))*1e-3); % [as/d]
             eopout_line{21}{curline+mjdind} = sprintf('%9.6f',yprat(mjdp==mjd(mjdind))*1e-3); % [as/d]     
             eopout_line{25}{curline+mjdind} = sprintf('%9.6f',xprat_e(mjdp==mjd(mjdind))*1e-3); % [as/d]
             eopout_line{26}{curline+mjdind} = sprintf('%9.6f',yprat_e(mjdp==mjd(mjdind))*1e-3); % [as/d]   
           else
             eopout_line{20}{curline+mjdind} = sprintf('%9s',filler);
             eopout_line{21}{curline+mjdind} = sprintf('%9s',filler); 
             eopout_line{25}{curline+mjdind} = sprintf('%9s',filler);
             eopout_line{26}{curline+mjdind} = sprintf('%9s',filler);            
           end
       else
           eopout_line{2}{curline+mjdind} = sprintf('%9s',filler);
           eopout_line{3}{curline+mjdind} = sprintf('%9s',filler);
           eopout_line{7}{curline+mjdind} = sprintf('%9s',filler);
           eopout_line{8}{curline+mjdind} = sprintf('%9s',filler);
           eopout_line{20}{curline+mjdind} = sprintf('%9s',filler);
           eopout_line{21}{curline+mjdind} = sprintf('%9s',filler); 
           eopout_line{25}{curline+mjdind} = sprintf('%9s',filler);
           eopout_line{26}{curline+mjdind} = sprintf('%9s',filler);
       end
       if indu(mjdind)
           if ~isnan(utfin(1))
               eopout_line{4}{curline+mjdind} = sprintf('%10.7f',utfin(mjdu==mjd(mjdind))*1e-3); % [s]
               eopout_line{9}{curline+mjdind} = sprintf('%10.7f',utfin_e(mjdu==mjd(mjdind))*1e-3); % [s]
           else
               eopout_line{4}{curline+mjdind} = sprintf('%10s',filler);
               eopout_line{9}{curline+mjdind} = sprintf('%10s',filler);
           end
           if ~isnan(lod(1))
             eopout_line{22}{curline+mjdind} = sprintf('%10.7f',lod(mjdu==mjd(mjdind))*1e-3); % [s/d]
             eopout_line{27}{curline+mjdind} = sprintf('%10.7f',utrat_e(mjdu==mjd(mjdind))*1e-3); % [s/d]   
           else
             eopout_line{22}{curline+mjdind} = sprintf('%10s',filler);
             eopout_line{27}{curline+mjdind} = sprintf('%10s',filler);     
           end
       else
           eopout_line{4}{curline+mjdind} = sprintf('%10s',filler);
           eopout_line{9}{curline+mjdind} = sprintf('%10s',filler);
           eopout_line{22}{curline+mjdind} = sprintf('%10s',filler);
           eopout_line{27}{curline+mjdind} = sprintf('%10s',filler);        
       end
       if indn(mjdind)
           if ~isnan(dXfin(1))
               eopout_line{5}{curline+mjdind} = sprintf('%6.3f',dXfin(mjdn==mjd(mjdind))); % [mas]
               eopout_line{6}{curline+mjdind} = sprintf('%6.3f',dYfin(mjdn==mjd(mjdind))); % [mas]
               eopout_line{10}{curline+mjdind} = sprintf('%6.3f',dXfin_e(mjdn==mjd(mjdind))); % [mas]
               eopout_line{11}{curline+mjdind} = sprintf('%6.3f',dYfin_e(mjdn==mjd(mjdind))); % [mas]
           else 
               eopout_line{5}{curline+mjdind} = sprintf('%6s',filler);
               eopout_line{6}{curline+mjdind} = sprintf('%6s',filler);
               eopout_line{10}{curline+mjdind} = sprintf('%6s',filler);
               eopout_line{11}{curline+mjdind} = sprintf('%6s',filler);
           end
           if ~isnan(dXrat(1))
               eopout_line{23}{curline+mjdind} = sprintf('%6.3f',dXrat(mjdn==mjd(mjdind))); % [mas/d]
               eopout_line{24}{curline+mjdind} = sprintf('%6.3f',dYrat(mjdn==mjd(mjdind))); % [mas/d]
               eopout_line{28}{curline+mjdind} = sprintf('%6.3f',dXrat_e(mjdn==mjd(mjdind))); % [mas/d]
               eopout_line{29}{curline+mjdind} = sprintf('%6.3f',dYrat_e(mjdn==mjd(mjdind))); % [mas/d] 
           else
               eopout_line{23}{curline+mjdind} = sprintf('%6s',filler);
               eopout_line{24}{curline+mjdind} = sprintf('%6s',filler);
               eopout_line{28}{curline+mjdind} = sprintf('%6s',filler);
               eopout_line{29}{curline+mjdind} = sprintf('%6s',filler);
           end
       else
           eopout_line{5}{curline+mjdind} = sprintf('%6s',filler);
           eopout_line{6}{curline+mjdind} = sprintf('%6s',filler);
           eopout_line{10}{curline+mjdind} = sprintf('%6s',filler);
           eopout_line{11}{curline+mjdind} = sprintf('%6s',filler);
           eopout_line{23}{curline+mjdind} = sprintf('%6s',filler);
           eopout_line{24}{curline+mjdind} = sprintf('%6s',filler);
           eopout_line{28}{curline+mjdind} = sprintf('%6s',filler);
           eopout_line{29}{curline+mjdind} = sprintf('%6s',filler);
       end
       eopout_line{12}{curline+mjdind} = sprintf('%7.2f',wrms_pfr);
       eopout_line{13}{curline+mjdind} = sprintf('%6s',filler);
       eopout_line{14}{curline+mjdind} = sprintf('%6s',filler);
       eopout_line{15}{curline+mjdind} = sprintf('%6s',filler);
       eopout_line{16}{curline+mjdind} = sprintf('%6s',filler);
       eopout_line{17}{curline+mjdind} = sprintf('%6.0f',num_obs);
       eopout_line{18}{curline+mjdind} = sprintf('%6s',IVSsesnam);
       eopout_line{19}{curline+mjdind} = sprintf('%5.2f',ses_dur);
       eopout_line{30}{curline+mjdind} = sprintf('%-64s',string(netID));
       
    end        
    curline = curline+mjdind;
end

% sort output cell array
[~,idx] = sort(eopout_line{1});
for n = 1:size(eopout_line,2)
     eopout_line{n} = eopout_line{n}(idx);
end

if flag_intensive
     eopout_name=[outfile,'.eopi'];
else
     eopout_name=[outfile,'.eoxy'];
end

mod_date = datestr(datetime());
fid = fopen(eopout_name,'w');
fprintf(fid, ['# Last updated: %s \n',...
                        '# Earth orientation parameters from the VLBI solution \n',...
                        '# VLBI VieVS solution \n',...
                        '# Analysis center: VIE (TU Wien, Austria)\n',...
                        '# \n',...
                        '# Time argument: UTC \n',...
                        '# Nutation angles are wrt IAU 2006 nutation/precession \n',...
                        '# \n',...
                        '#    EOXY file contains the series of estimates of the Earth orientation \n',...
                        '# parameters obtained from processing VLBI experiments. \n',...
                        '# \n',...
                        '# Format: \n',...
                        '# \n',...
                        '# The line which starts from # is considered as a comment \n',...
                        '# \n',...
                        '# Field Columns Format Units   Meaning \n',...
                        '# \n',...
                        '#   1   1-13    F12.6  days    Modified Julian date of the UTC time tag \n',...
                        '#   2   14-23   F9.6   arcsec  The estimate of X pole coordinate \n',...
                        '#   3   24-33   F9.6   arcsec  The estimate of Y pole coordinate \n',...
                        '#   4   34-44   F10.7  sec     The UT1-UTC function \n',...
                        '#   5   45-51   F6.3   mas     Celestial Pole Offset dX with \n',...
                        '#                                 respect to IAU 2006 nutation/precession \n',...
                        '#   6   52-58   F6.3   mas     Celestial Pole Offset dY with \n',...
                        '#                                 respect to IAU 2006 nutation/precession \n',...
                        '#   7   59-68   F9.6   arcsec  Formal uncertainty of X pole coordinate \n',...
                        '#   8   69-78   F9.6   arcsec  Formal uncertainty of Y pole coordinate \n',...
                        '#   9   79-89   F10.7  sec     Formal uncertainty of UT1-UTC function \n',...
                        '#  10   90-96   F6.3   mas     Formal uncertainty of nutation dX \n',...
                        '#  11   97-103  F6.3   mas     Formal uncertainty of nutation dY \n',...
                        '#  12   104-111 F7.2   psec    Weighted root mean square of postfit residual of \n',...
                        '#                                 the solution \n',...
                        '#  13   112-118 F6.4   --      Correlation between the estimates of X-pole \n',...
                        '#                                 positions and Y-pole position \n',...
                        '#  14   119-125 F6.4   --      Correlation between the estimates of X-pole \n',...
                        '#                                 positions and UT1-TAI angle \n',...
                        '#  15   126-132 F6.4   --      Correlation between the estimates of Y-pole \n',...
                        '#                                 positions and UT1-TAI angle \n',...
                        '#  16   133-139 F6.4   --      Correlation between the estimates of nutation \n',...
                        '#                                 dX and nutation dY \n',...
                        '#  17   140-146 I6     --      Number of used observations in the session \n',...
                        '#  18   147-153 A6     --      IVS session code \n',...
                        '#  19   154-159 F5.2   hours   Session duration \n',...
                        '#  20   160-169 F9.6   asc/day Rate of X pole coordinate \n',...
                        '#  21   170-179 F9.6   asc/day Rate of Y pole coordinate \n',...
                        '#  22   180-190 F10.7  sec     Length of day \n',...
                        '#  23   191-197 F6.3   mas/day Rate of dX \n',...
                        '#  24   198-204 F6.3   mas/day Rate of dY \n',...
                        '#  25   205-214 F9.6   asc/day Formal uncertainty of X pole coordinate rate \n',...
                        '#  26   215-224 F9.6   asc/day Formal uncertainty of Y pole coordinate rate \n',...
                        '#  27   225-235 F10.7  sec     Formal uncertainty of length of day \n',...
                        '#  28   236-242 F6.3   mas/day Formal uncertainty of dX rate \n',...
                        '#  29   243-249 F6.3   mas/day Formal uncertainty of dY rate \n',...
                        '#  30   250-315 A64    --      Network ID. The alphabetically ordered sequence \n',...
                        '#                              of two-letter IVS station identifiers. Only those \n',...
                        '#                              stations which provided observations used in the \n',...
                        '#                              solution are listed. The station names are defined \n',...
                        '#                              in the IVS document ivscontrol/ns-codes.txt \n',...
                        '# \n',...
                        '# If a given parameter was not estimated a filler, -0, is placed. \n',...
                        '# \n'],mod_date);

    for n = 1:length(eopout_line{1})
          fprintf(fid,[repmat('%s ',1,30),'\n'],string(cellfun(@(x) x(n), eopout_line))');
    end
    fclose(fid);
