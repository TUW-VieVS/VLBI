% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% This is a tool to generate a textfile with the EOP results from the
% vie_lsm output. It should be executed in the OUT directory. The 
% function reads the actual process_list from the WORK directory and 
% creates one textfile (eop_subdir.txt) for all sessions on that list.
% The x_ files in DATA/LEVEL3 are used.
% If nutation and erp were estimated for different times, linear
% interpolation is done.
% The total values (colums 2-6) do not contain the model values for high 
% frequency variations.
% tot = a priori + estimated
% coded by Lucia, 08/2009
% adopted for VieVS_V1c by Lucia, 08/2010
% modified by Sigrid Boehm 09/2011
% new tidal variations in dUT1 before/After interpolation by Lucia, 04/2013
% variable output filename possible -> give varargin{1} full path+filename.
% 	by Matthias, 03/2014
% process_list can also be a saved process list (.mat file)
% 	by Matthias, 03/2014
% modification so that the program can be used to output ut1 from
% Intensives
%   by Johannes, 24/03/2014
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function eop_out_int(process_list,subdir,varargin)

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

if nargin>2
	outfile=varargin{1};
else
	outfile=['../OUT/eop_',subdir,'.txt'];
end

fid = fopen(outfile,'wt');
fprintf(fid,'%%************************************************************\n');
fprintf(fid,'%% Columns:\n%%\t 1     .... mjd\n%%\t 2-6   .... total values(x,y,ut,dX,dY)\n%%\t 7-11  .... a priori EOP (input in vie_mod)\n%%\t 11-16 .... estimated values\n%%\t 17-21 .... error of estimation\n%%\t 22-24 .... high frequency (subdaily) ERP corrections\n%%\n'); 
fprintf(fid,'%% all units in mas resp. ms (dut1)\n');
fprintf(fid,'%%************************************************************\n');
fprintf(fid,'%%    MJD       xpol        ypol        dut1        dX          dY          x_apr       y_apr       ut_apr      dX_apr      dY_apr      x_est       y_est       dut1_est    dX_est      dY_est      x_err       y_err       ut_err      dX_err      dY_err      x_hf        y_hf        ut_hf\n');
fprintf(fid,'%%\n');

for j = 1:pl(1)
	sname = process_list(j,end-13:end);

	if ~(exist(strcat('../DATA/LEVEL3/',subdir,'/x_',num2str(sname),'.mat'), 'file'))
	   clearvars -except subdir pl process_list fid j;
	   continue
	end
	load(strcat('../DATA/LEVEL3/',subdir,'/x_',num2str(sname)));
	load(strcat('../DATA/LEVEL3/',subdir,'/opt_',num2str(sname)));
	load(strcat('../DATA/LEVEL1/',subdir,'/',num2str(sname),'_parameter'));

	if opt_.xpol.model==1
		mjdp = x_.xpol.mjd;
	else
		mjdp = [];
	end

	if (length([x_.dut1.mjd]) ~= length(mjdp))
%		mjdu = x_.dut1.mjd; % Johannes
        mjdu = (opt_.first_scan + opt_.last_scan)/2;
	else
		mjdu = mjdp;
	end

	% nutation if estimated
	if opt_.nutdx.model==1
	   mjdn   = x_.nutdx.mjd;
	else
	   mjdn   = [];
	end

	% make one vector with all mjd values
	mjd = [mjdp mjdu mjdn];
	mjd = sort(mjd);
	mjd(diff(mjd)==0)=[];
	mjd(diff(mjd)==0)=[];

	% find the indices of the estimated intervalls in mjd
	for i=1:length(mjdp)
		indp(i) = find(mjd==mjdp(i));
	end
	for i=1:length(mjdu)
		indu(i) = find(mjd==mjdu(i));
	end
	for i=1:length(mjdn)
		indn(i) = find(mjd==mjdn(i));
	end

	% create NaN vectors for each component with the length of mjd
	% and fill in the values at the corresponding mjd place
	xp   = NaN(1,length(mjd));    
	xp_e = NaN(1,length(mjd));    
	yp   = NaN(1,length(mjd));       
	yp_e = NaN(1,length(mjd));     
	ut   = NaN(1,length(mjd));     
	ut_e = NaN(1,length(mjd));     
	dX   = NaN(1,length(mjd));      
	dX_e = NaN(1,length(mjd));     
	dY   = NaN(1,length(mjd));      
	dY_e = NaN(1,length(mjd)); 
	if opt_.xpol.model==1
	   xp(indp) = x_.xpol.val;  
	   xp_e(indp) = x_.xpol.mx;
	   yp(indp) = x_.ypol.val;
	   yp_e(indp) = x_.ypol.mx; 
	end
	if opt_.dut1.model==1
	   ut(indu) = x_.dut1.val(1);   % Johannes
	   ut_e(indu) = x_.dut1.mx(1);  % Johannes
	end
	if opt_.nutdx.model==1
	   dX(indn) = x_.nutdx.val; 
	   dX_e(indn) = x_.nutdx.mx; 
	   dY(indn) = x_.nutdy.val;  
	   dY_e(indn) = x_.nutdy.mx;   
	end 

	% high frequency corrections
	LEAP  = tai_utc(mjd);                   % get difference TAI-UTC [s]
	TT = mjd + (32.184 + LEAP')./86400;     % [MJD TT]
	% set parameters for triaxiality models to zero to get only the ocean tidal
	% influence
	% parameter.vie_mod.lib_pm = 0;
	% parameter.vie_mod.lib_ut = 0;
	[xhf,yhf,uthf] = eophf(TT',parameter);   % [rad, sec]
	xhf = rad2mas(xhf); yhf = rad2mas(yhf); uthf = uthf *1e3; %[mas, ms]

	% get a priori values
	MJDeop =  parameter.eop.mjd;
	XPeop  = (parameter.eop.xp)*1e3;  % [mas]
	YPeop  = (parameter.eop.yp)*1e3;  % [mas]
	UT1eop = (parameter.eop.ut1)*1e3; % [ms]
	dXeop  = (parameter.eop.dX)*1e3;  % [mas]
	dYeop  = (parameter.eop.dY)*1e3;  % [mas]   

    % subtraction of tidal variations (Defraigne and Smits) in dUT1 before
    % interpolation
    if parameter.vie_mod.tidalUT == 1
        disp('remove tidal UT')
        taiut    = tai_utc(MJDeop); 
        MJDTTeop = MJDeop + (32.184 + taiut)/86400;
        %         UT1corr  = tver2000(MJDTTeop);  % [sec]
        if parameter.vie_mod.tidalUT35 ==1
            par35=1;
        else
            par35=2;
        end
        UT1corr  = rg_zont2(MJDTTeop,par35);  % [sec]
        UT1eop   = UT1eop - UT1corr*1e3;  % [ms]
    end 

	if strcmp(parameter.eop.interp,'linear')
		DUT1  = interp1(MJDeop,UT1eop,mjd,'linear','extrap');  
		XP    = interp1(MJDeop, XPeop,mjd,'linear','extrap');  
		YP    = interp1(MJDeop, YPeop,mjd,'linear','extrap');  
		DXap  = interp1(MJDeop, dXeop,mjd,'linear','extrap');
		DYap  = interp1(MJDeop, dYeop,mjd,'linear','extrap');
	else % lagrange interpolation
		% Subtraction of tidal variations (Defraigne and Smits) in dut1
		DUT1  = (lagint4v(MJDeop,UT1eop,mjd))';  
		XP    = (lagint4v(MJDeop, XPeop,mjd))';  
		YP    = (lagint4v(MJDeop, YPeop,mjd))';  
		DXap  = (lagint4v(MJDeop, dXeop,mjd))';
		DYap  = (lagint4v(MJDeop, dYeop,mjd))';
	end

    % re-add tidal variation in dUT1
    if parameter.vie_mod.tidalUT == 1
        disp('re-add tidal UT')
        %corrUT1 = tver2000(TT);     % [sec]
        corrUT1 = rg_zont2(TT,par35);  % [sec]
        DUT1    = DUT1 + corrUT1'*1e3;   % [ms]
    end

	if parameter.vie_mod.eophf ==1;
		inclhf = 'yes';
	else
		inclhf = 'no';
	end

	if parameter.vie_mod.dXdY == 1
	   inclanut = 'yes';
	else
	   inclanut = 'no';
	end


	% calculate total EOP values (tot = estimated + a priori)
	xptot = xp + XP;   % [mas]
	yptot = yp + YP;   % [mas]
	uttot = ut + DUT1; % [ms]
	dXtot = dX + DXap; % [mas]
	dYtot = dY + DYap; % [mas]


	% create file
	% unit: mas/ms
	out = [mjd;xptot;yptot;uttot;dXtot;dYtot;XP;YP;DUT1;DXap;DYap;xp;yp;ut;dX;dY;xp_e;yp_e;ut_e;dX_e;dY_e;xhf';yhf';uthf'];

	% fid = fopen(strcat('../OUT/eop_',num2str(sname),'.txt'),'wt');
	% fprintf(fid,'%%************************************************************\n');
	% fprintf(fid,strcat('%% EOP out table of lsm file x_',num2str(sname),'\n'));
	% fprintf(fid,'%% Columns:\n%%\t 1     .... mjd\n%%\t 2-6   .... total values(x,y,ut,dX,dY)\n%%\t 7-11  .... a priori EOP (input in vie_mod)\n%%\t 12-16 .... estimated values\n%%\t 17-21 .... error of estimation\n%%\t 22-24 .... high frequency (subdaily) ERP corrections\n%%\n'); 
	% fprintf(fid,'%% all units in mas resp. ms (dut1)\n');
	% fprintf(fid,strcat('%% subdaily high frequency variations included in the estimations: ',inclhf,'\n')); 
	% fprintf(fid,strcat('%% a priori nutation corrections were applied:',inclanut,'\n'));
	% fprintf(fid,'%%************************************************************\n');
	% fprintf(fid,'%%    MJD       xpol        ypol        dut1        dX          dY          x_apr       y_apr       ut_apr      dX_apr      dY_apr      x_est       y_est       dut1_est    dX_est      dY_est      x_err       y_err       ut_err      dX_err      dY_err      x_hf        y_hf        ut_hf\n');
	% fprintf(fid,'%%\n');

	fprintf(fid,' %5.4f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f  ',out);
	fprintf(fid,' %s\n',sname);
    % fclose(fid);
	clearvars -except subdir pl process_list fid j;

end
fclose(fid);