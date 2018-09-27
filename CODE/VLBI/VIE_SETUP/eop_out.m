% #########################################################################
% #     eop_out
% #########################################################################
%
% DESCRIPTION
%   This is a tool to generate a textfile with the EOP results from the
%   vie_lsm output. It should be executed in the OUT directory. The 
%   function reads the actual process_list from the WORK directory and 
%   creates one textfile (eop_subdir_format.txt) for all sessions on that list.
% 
%   The x_ files in DATA/LEVEL3 are used.
%
%   If nutation and ERP were estimated for different times, common mjd vector is generated
%		and values estimated less frequently are filled out with NAN
%	for example, ERP hourly and nutation daily. 
%	In output you have NAN values for nutation on hourly intervals and 
%	  you have estimated nutation values at midnights only
%
%   The total values (colums 2-6) do not contain the model values for high 
%   frequency variations.
%   tot = a priori + estimated
%
%	Note! The set of parameters for triaxiality models to zero 
% 			to get the ocean tidal influence only
%
% UPDATES:
%	internal vector zeit should be updated when a new leap second is introduced
%
% CREATED  
%   2009/08     Lucia Plank
%
% REFERENCES
%
%
% COUPLING
% - azel2xyew	???
%
%   External calls:
%		/VIE_MOD/tai_utc.m
%		/VIE_MOD/eophf.m
%		/VIE_MOD/rad2mas.m
%		/VIE_MOD/rg_zont2.m
%		    /OUT/get_out_val.m
%			/OUT/eop_get_approx.m
%		    /OUT/weightedMeanFunction.m (in get_out_val.m)
%
% INPUT
% - process_list            : Process list containing all sessions which results should be written to the output file(s) 
% - subdir                  : Sub-directory in VieVS file structure
% - outfile (optional)      : Output filename and directory
% - output_mode (optional)  : Parameter to select between different output file formats, 
%								format: [flag_file_version_1, flag_file_version_2, flag_file_version_3]
%
% OUTPUT
%	creates txt-files with respect to input option
%	I. Regular data
%   There are three different options for the output format of the EOP file:
%       1.) output_mode(1) = 1: detailed EOP file (as before)
%       2.) output_mode(2) = 1: Sorted EOP file (sorted by time/date, multiple entries are thrown out)
%       3.) output_mode(3) = 1: VieVS specific format (can be used as a priori model in VieVS, 
%								5 values to edges are added using data form LEVEL1/parameter-structure)
%	II. Intensive session
%		output_mode(1) = 1: detailed EOP file (as before)
%
% CHANGES:
% - 2010-08-??, L. Plank: adopted for VieVS_V1c
% - 2011-09-??, S. Böhm: modified
% - 2013-04-??, L. Plank: new tidal variations in dUT1 before/After interpolation
% - 2014-03-??, M. Madzak: variable output filename possible -> give varargin{1} full path+filename
% - 2014-03-??, M. Madzak: process_list can also be a saved process list (.mat file)
% - 2015-08-19, A. Girdiuk: Possibility added to choose between different output options:
%                            1.) output_mode(1) = 1: detailed EOP file (as before)
%                            2.) output_mode(2) = 1: Sorted EOP file (sorted by time/date, multiple entries are thrown out)
%                            3.) output_mode(3) = 1: VieVS specific format (can be used as a priori model in VieVS,
%													 exctra-EOP values are added from finals)
% - 2015-09-14, A. Girdiuk: bug-fix
% - 2016-04-25, A. Girdiuk: bug-fix
% - 2016-08-08, A. Girdiuk, L. Plank : transformation to TT and leap-second case solved
% - 2016-09-19, A. Girdiuk: suited to make intensive sessions output
% - 2016-11-29, A. Girdiuk: suited to use a priori data from
%                           LEVEL1/parameter-structure for VieVS specific format (output_mode(3))
%							'The set of parameters for triaxiality models to zero' is uncommented
%							function description is updated
%                           errors are added to format (output_mode(2))
% - 2016-11-29, A. Girdiuk: re-viewed; bug-fix
% - 2018-03-13, A. Hellerschmied: function now supports arbitrary session names in process lists.
% - 2018-09-25, S. Boehm: if default tight constraints were used for EOP only one EOP value per session is written,
%                          referred to the middle of the session.

function eop_out(process_list, subdir, varargin)

% init.:
output_mode = [1,0,0]; % Default: Detailed EOP outout.
flag_intensive = 0;

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

if nargin == 2
    outfile=['../OUT/eop_',subdir]; % default
    
elseif nargin == 3 % outfile name and directory is defined in the function input
    outfile = varargin{1}(1:(length(varargin{1})-4));
    
elseif nargin == 4 % An output option is defined in the function input
    outfile     = varargin{1}(1:(length(varargin{1})-4));
    output_mode = logical(varargin{2});
elseif nargin == 5
    outfile     = varargin{1}(1:(length(varargin{1})-4));
    flag_intensive= logical(varargin{3});
else
    % ERROR!
    fprintf(1, 'ERROR in eop_out.m: Wrong number of input arguments!');
end

%output_mode = [0,0,0]
%   	1	   [1,0,0]	    detailed (standard)
%	    2	   [0,1,0]	    sorted
%	   	3	   [0,0,1]	 	vievs format (sorted + filled)

if flag_intensive
	outfile_int=[outfile,'_intensive','.txt'];
	fid_int = fopen(outfile_int,'wt');
	fprintf(fid_int,'%%************************************************************\n');
	fprintf(fid_int,'%% Columns:\n%%\t 1     .... mjd\n%%\t 2     .... total values\n%%\t 3     .... a priori EOP (input in vie_mod)\n%%\t 4     .... estimated values\n%%\t 5     .... error of estimation\n%%\t 6     .... high frequency (subdaily) ERP corrections\n%%\n'); 
	fprintf(fid_int,'%% all units in mas resp. ms (dut1)\n');
	fprintf(fid_int,'%%************************************************************\n');
	fprintf(fid_int,'%%    MJD        dut1        ut_apr      dut1_est     ut_err      ut_hf\n');
	fprintf(fid_int,'%%\n');
else
    % 1
    if output_mode(1)
        outfile_1=[outfile,'_detailed','.txt'];
        fid_1 = fopen(outfile_1,'wt');
        fprintf(fid_1,'%%************************************************************\n');
        fprintf(fid_1,'%% Columns:\n%%\t 1     .... mjd\n%%\t 2-6   .... total values(x,y,ut,dX,dY)\n%%\t 7-11  .... a priori EOP (input in vie_mod)\n%%\t 12-16 .... estimated values\n%%\t 17-21 .... error of estimation\n%%\t 22-24 .... high frequency (subdaily) ERP corrections\n%%\n'); 
        fprintf(fid_1,'%% all units in mas resp. ms (dut1)\n');
        fprintf(fid_1,'%%************************************************************\n');
        fprintf(fid_1,'%%    MJD       xpol        ypol        dut1        dX          dY          x_apr       y_apr       ut_apr      dX_apr      dY_apr      x_est       y_est       dut1_est    dX_est      dY_est      x_err       y_err       ut_err      dX_err      dY_err      x_hf        y_hf        ut_hf\n');
        fprintf(fid_1,'%%\n');
    end
    % 2
    if output_mode(2)
        outfile_2=[outfile,'_ordered','.txt'];
        fid_2 = fopen(outfile_2,'wt');
        fprintf(fid_2,'%%************************************************************\n');
        fprintf(fid_2,'%% Columns:\n%%\t 1     .... mjd\n%%\t 2-6   .... total values(x,y,ut,dX,dY)\n');
        fprintf(fid_2,'%%	 7-12  .... errors of estimation(x,y,ut,dX,dY)\n');
        fprintf(fid_2,'%% all units in mas resp. ms (dut1)\n');
        fprintf(fid_2,'%%************************************************************\n');
        fprintf(fid_2,'%%    MJD       xpol      ypol     dut1       dX     dY         xp_err    yp_err  dut1_err  dX_err dY_err\n');
    end
    % 3
    if output_mode(3)
        outfile_3=[outfile,'_vievs-format','.txt'];
        fid_3 = fopen(outfile_3,'wt');
        fprintf(fid_3,'%%************************************************************\n');
        fprintf(fid_3,'%% Columns:\n%%\t 1     .... mjd\n%%\t 2-6   .... total values(x,y,ut,dX,dY)\n');
        fprintf(fid_3,'%% all units in mas resp. ms (dut1)\n');
        fprintf(fid_3,'%%************************************************************\n');
        fprintf(fid_3,'%% MJD  xpol  ypol  dut1  dX  dY\n');
    end
end

%%% initialization
mjd_keep  = [];
xp_keep   = []; xp_e_keep   = [];  X_keep    = [];  X_e_keep = []; 
yp_keep   = []; yp_e_keep   = [];  Y_keep    = [];  Y_e_keep = []; 
Dut1_keep = []; Dut1_e_keep = [];
%%%

box=[];
for j = 1:pl(1)
    
    % Get session name:
    % Unix:
    ind_unix = strfind(process_list(j,:), '/');
    % DOS:
    ind_dos = strfind(process_list(j,:), '\');
    ind = max([ind_unix, ind_dos]);
    sname = process_list(j,ind+1:end);

	load(strcat('../DATA/LEVEL3/',subdir,'/x_',sname));
	load(strcat('../DATA/LEVEL3/',subdir,'/opt_',sname));
	load(strcat('../DATA/LEVEL1/',opt_.level1OutDir,'/',sname,'_parameter'));

	if parameter.lsmopt.xpol.model==1
        if opt_.xpol.coef > 1.0e-4
		   mjdp = x_.xpol.mjd;
        else
           mjdp = (opt_.first_scan + opt_.last_scan)/2; 
        end
	else
		mjdp = [];
	end

	if parameter.lsmopt.dut1.model==1
        if flag_intensive || (opt_.dut1.coef <= 1.0e-4)
            mjdu = (opt_.first_scan + opt_.last_scan)/2;
        else
            mjdu = x_.dut1.mjd;
        end
	else
		mjdu = [];
	end

	% nutation if estimated
	if parameter.lsmopt.nutdx.model==1
       if opt_.nutdx.coef > 1.0e-4
	      mjdn   = x_.nutdx.mjd;
       else
          mjdn = (opt_.first_scan + opt_.last_scan)/2; 
       end
	else
	   mjdn   = [];
	end

	% make one vector with all mjd values
	mjd = [mjdp mjdu mjdn];
	mjd = sort(mjd);
	mjd(diff(mjd)==0)=[];
    
	%initialization
    indp = ismember(mjd , mjdp);
    indu = ismember(mjd , mjdu);
    indn = ismember(mjd , mjdn);
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
    
	if parameter.lsmopt.xpol.model==1
	   xp(indp)   = x_.xpol.val(1:length(mjdp));  
	   xp_e(indp) = x_.xpol.mx(1:length(mjdp));
	   yp(indp)   = x_.ypol.val(1:length(mjdp));
	   yp_e(indp) = x_.ypol.mx(1:length(mjdp)); 
	end
	if parameter.lsmopt.dut1.model==1
	   ut(indu)   = x_.dut1.val(1:length(mjdu));   
	   ut_e(indu) = x_.dut1.mx(1:length(mjdu)); 
	end
	if parameter.lsmopt.nutdx.model==1
	   dX(indn)   = x_.nutdx.val(1:length(mjdn)); 
	   dX_e(indn) = x_.nutdx.mx(1:length(mjdn)); 
	   dY(indn)   = x_.nutdy.val(1:length(mjdn));  
	   dY_e(indn) = x_.nutdy.mx(1:length(mjdn));   
	end 

	% high frequency corrections
	[LEAP,zeit]  = tai_utc(mjd,1);          % get difference TAI-UTC [s]
	TT = mjd + (32.184 + LEAP')./86400;     % [MJD TT]
	% set parameters for triaxiality models to zero to get only 
	% the ocean tidal influence
	parameter.vie_mod.lib_pm = 0;
	parameter.vie_mod.lib_ut = 0;
	[xhf,yhf,uthf] = eophf(TT',parameter);   % [rad, sec]
	xhf = rad2mas(xhf); yhf = rad2mas(yhf); uthf = uthf *1e3; %[mas, ms]

    [XP,YP,DUT1,DXap,DYap,MJDeop] = eop_get_approx(parameter,mjd,TT);
    
    
    box(j,1) = min(MJDeop);
    box(j,2) = max(MJDeop);

    % calculate total EOP values (tot = estimated + a priori)

    uttot = ut + DUT1; % [ms]

	if ~isempty(intersect(zeit,mjd)) && any(diff(LEAP)~=0)
        [~,~,IB]=intersect(zeit,mjd);
    	uttot(IB:end) = ut(IB:end) + DUT1(IB:end) +1e3; % [ms]
    end

    xptot = xp + XP;   % [mas]
    yptot = yp + YP;   % [mas]
        %
    dXtot = dX + DXap; % [mas]
	dYtot = dY + DYap; % [mas]

    % Int. output
    if flag_intensive
        out = [mjd;uttot(1);DUT1;ut(1);ut_e(1);uthf'];
		fprintf(fid_int,' %5.4f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f\n',out);
    end
    
	if output_mode(1)
		%%%	> 1.	The detailed EOP data list
		% create file
		% unit: mas/ms
		out = [mjd;xptot;yptot;uttot;dXtot;dYtot;XP;YP;DUT1;DXap;DYap;xp;yp;ut;dX;dY;xp_e;yp_e;ut_e;dX_e;dY_e;xhf';yhf';uthf'];
		fprintf(fid_1,' %5.4f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f\n',out);
		%%% <.1
	end
%%%
	if output_mode(2) || output_mode(3)
		%%% > 2. The ordered EOP data list
		%%%

		mjd_keep  = [mjd_keep;   mjd'];
		xp_keep   = [ xp_keep; xptot'];
		yp_keep   = [ yp_keep; yptot'];
   		Dut1_keep = [Dut1_keep;uttot'];
  		X_keep    = [X_keep;   dXtot'];
   		Y_keep    = [Y_keep;   dYtot'];
   		% errors
    	xp_e_keep   = [ xp_e_keep; xp_e'];
   		yp_e_keep   = [ yp_e_keep; yp_e'];
   		Dut1_e_keep = [Dut1_e_keep;ut_e'];
		X_e_keep    = [X_e_keep;   dX_e'];
		Y_e_keep    = [Y_e_keep;   dY_e'];
	end

end

if output_mode(2) || output_mode(3)

	outMjd = sort(unique(mjd_keep));
	outxp  = zeros(length(outMjd),1); % out = sorted, one per mjd
	outyp  = outxp;
	outDut1= outxp;
	outX   = outxp;
	outY   = outxp;

	for iMjd=1:length(outMjd)
		foundLines = mjd_keep==outMjd(iMjd);
		if sum(foundLines)==1
			outxp(iMjd)  =   xp_keep(foundLines);   outxp_err(iMjd)  =   xp_e_keep(foundLines);
			outyp(iMjd)  =   yp_keep(foundLines);   outyp_err(iMjd)  =   yp_e_keep(foundLines);
			outDut1(iMjd)= Dut1_keep(foundLines);   outDut1_err(iMjd)= Dut1_e_keep(foundLines);
			outX(iMjd)   =    X_keep(foundLines);   outX_err(iMjd)   =    X_e_keep(foundLines);
			outY(iMjd)   =    Y_keep(foundLines);   outY_err(iMjd)   =    Y_e_keep(foundLines);
        else
           [outxp(iMjd),outxp_err(iMjd)]     = get_out_val(xp_keep(foundLines),    xp_e_keep(foundLines));
           [outyp(iMjd),outyp_err(iMjd)]     = get_out_val(yp_keep(foundLines),    yp_e_keep(foundLines));
           [outDut1(iMjd),outDut1_err(iMjd)] = get_out_val(Dut1_keep(foundLines),Dut1_e_keep(foundLines));
           [outX(iMjd),outX_err(iMjd)]       = get_out_val(X_keep(foundLines),      X_e_keep(foundLines));
           [outY(iMjd),outY_err(iMjd)]       = get_out_val(Y_keep(foundLines),      Y_e_keep(foundLines));
		end
	end

	if output_mode(2)
	%%%	> 2.	The detailed EOP data list
		out_2 = [outMjd';outxp';outyp';outDut1';outX';outY';...
                          outxp_err;outyp_err;outDut1_err;outX_err;outY_err];
		fprintf(fid_2,'%7.2f %9.3f %9.3f %10.4f %7.4f %7.4f  %9.3f %9.3f %10.4f %7.4f %7.4f\n',out_2);
	%%% <2.
	end

	if output_mode(3)
	%%% > 3. VieVS specific format
        [~,ind_min]=min(box(:,1));% = min(MJDeop);
        [~,ind_max]=max(box(:,2));% = max(MJDeop);
        sname = process_list(ind_min,6:end);
        load(strcat('../DATA/LEVEL1/',opt_.level1OutDir,'/',num2str(sname),'_parameter'));
        % get a priori values
        MJDeop =  parameter.eop.mjd;
        XPeop  = (parameter.eop.xp)*1e3;  % [mas]
        YPeop  = (parameter.eop.yp)*1e3;  % [mas]
        UT1eop = (parameter.eop.ut1)*1e3; % [ms]
        dXeop  = (parameter.eop.dX)*1e3;  % [mas]
        dYeop  = (parameter.eop.dY)*1e3;  % [mas]  
        
        if find(MJDeop==min(outMjd))<5
            disp('Warning: number edge points less than 5: beginning of time span\n')
        end
        
        MJDeop_add_edge_min = MJDeop(1:find(MJDeop==min(outMjd))-1);
        XPeop_add_edge_min  = XPeop(1:find(MJDeop==min(outMjd))-1);
        YPeop_add_edge_min  = YPeop(1:find(MJDeop==min(outMjd))-1);
        UT1eop_add_edge_min = UT1eop(1:find(MJDeop==min(outMjd))-1);
        dXeop_add_edge_min  = dXeop(1:find(MJDeop==min(outMjd))-1);
        dYeop_add_edge_min  = dYeop(1:find(MJDeop==min(outMjd))-1);

        sname = process_list(ind_max,6:end);
        load(strcat('../DATA/LEVEL1/',opt_.level1OutDir,'/',num2str(sname),'_parameter'));
        % get a priori values
        MJDeop =  parameter.eop.mjd;
        XPeop  = (parameter.eop.xp)*1e3;  % [mas]
        YPeop  = (parameter.eop.yp)*1e3;  % [mas]
        UT1eop = (parameter.eop.ut1)*1e3; % [ms]
        dXeop  = (parameter.eop.dX)*1e3;  % [mas]
        dYeop  = (parameter.eop.dY)*1e3;  % [mas] 
        
        if find(MJDeop==max(outMjd))<5
            disp('Warning: number edge points less than 5: end of time span\n')
        end
        
        MJDeop_add_edge_max = MJDeop(find(MJDeop==max(outMjd))+1:end);
        XPeop_add_edge_max  = XPeop(find(MJDeop==max(outMjd))+1:end);
        YPeop_add_edge_max  = YPeop(find(MJDeop==max(outMjd))+1:end);
        UT1eop_add_edge_max = UT1eop(find(MJDeop==max(outMjd))+1:end);
        dXeop_add_edge_max  = dXeop(find(MJDeop==max(outMjd))+1:end);
        dYeop_add_edge_max  = dYeop(find(MJDeop==max(outMjd))+1:end);
        
        % pooling edges and main set
        aEOP_mjd  = [MJDeop_add_edge_min;outMjd; MJDeop_add_edge_max];
        aEOP_xp   = [XPeop_add_edge_min; outxp;  XPeop_add_edge_max];
        aEOP_yp   = [YPeop_add_edge_min; outyp;  YPeop_add_edge_max];
        aEOP_dut1 = [UT1eop_add_edge_min;outDut1;UT1eop_add_edge_max];
        aEOP_X    = [dXeop_add_edge_min; outX;   dXeop_add_edge_max];
        aEOP_Y    = [dYeop_add_edge_min; outY;   dYeop_add_edge_max];
        
         out_3 = [aEOP_mjd';aEOP_xp';aEOP_yp';aEOP_dut1';aEOP_X';aEOP_Y'];
 
         fprintf(fid_3,'%7.2f %9.3f %9.3f %10.4f %7.4f %7.4f\n',out_3);
        
	end

end
%
%	close existing files
%
if output_mode(1)
	fclose(fid_1);
end
%%% <.1
if output_mode(2)
	fclose(fid_2);	
end
%%% <.2
if output_mode(3)
	fclose(fid_3);
end
%%% <.3

if flag_intensive
    fclose(fid_int);
end
