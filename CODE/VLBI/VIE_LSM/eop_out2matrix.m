function eop_struct=eop_out2matrix(x_, opt_, parameter, maxMjdEop)
% adopted from eop_out.m

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% This is a tool to generate a textfile with the EOP results from the
% vie_lsm output. It should be executed in the OUT directory. The 
% function reads the actual process_list from the WORK directory and 
% creates one textfile (eop_sessionname.txt) for each session on that list.
% The x_ files in DATA/LEVEL3 are used.
% If nutation and erp were estimated for different times, linear
% interpolation is done.
% The total values (colums 2-6) are the instantaneous values, including  
% high frequency variations.
%   tot = a priori + estimated

% coded by Lucia, 08/2009
% adopted for VieVS_V1c by Lucia, 08/2010
% changed:  to function (output struct, no textfile)
%           by Matthias, 08/2010
% changed:  all variables are now interpolated independently to the maximum
%           number of mjds. If there are more than one variable with the
%           same max number of estimated values, the first one
%           (xpol->ypol->dut1->nutx->nuty) is used.
%           by Matthias 11/2010
% changed:  When no estimation was done in VieVS, this function can be 
%           called with x_=[] and opt_=[]: Then only a-priori values are
%           calculated correctly, estimates and errors are set to 0.
%           by Matthias 01/2011
% changed:  DUT1corr*1e3
%           by Lucia, 06/2011
% new tidal variations in dUT1 before/After interpolation by Lucia, 04/2013
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% define epochs at which all eops are to be interpolated (min interval)
mjd=maxMjdEop;

% if both x_ and opt_ exist (and are not empty) := normal with estimation:
if ~isempty(x_) && ~isempty(opt_)
    
    if opt_.xpol.model==1
        mjde = x_.xpol.mjd;
    elseif opt_.ypol.model==1
        mjde = x_.ypol.mjd;
    elseif opt_.dut1.model==1
        mjde = x_.dut1.mjd;
    elseif opt_.nutdx.model==1
        mjde = x_.nutdx.mjd;
    elseif opt_.nutdy.model==1
        mjde = x_.nutdy.mjd;
    else
        error('error writing eop data')
    end
    
    
    if opt_.xpol.model==1
       mjdxpol=x_.xpol.mjd;
       %mjdn   = x_.nutdx.mjd;
       xp   = x_.xpol.val;
       xp_e = x_.xpol.mx;
    else
       mjdxpol= [mjde(1);mjde(end)];
       xp     = [0; 0];
       xp_e   = [0; 0];
    end
    
    % ypol if estimated
    if opt_.ypol.model==1
       mjdypol=x_.ypol.mjd;
       %mjdn   = x_.nutdx.mjd;
       yp   = x_.ypol.val;
       yp_e = x_.ypol.mx;
    else
       mjdypol= [mjde(1);mjde(end)];
       yp     = [0; 0];
       yp_e   = [0; 0];
    end

    % dut1 if estimated
    if opt_.dut1.model==1
       mjddut1=x_.dut1.mjd;
       %mjdn   = x_.nutdx.mjd;
       ut   = x_.dut1.val;
       ut_e = x_.dut1.mx;
    else
       mjddut1= [mjde(1);mjde(end)];
       ut     = [0; 0];
       ut_e   = [0; 0];
    end

    % nutation dX if estimated
    if opt_.nutdx.model==1
       mjdnutdx=x_.nutdx.mjd;
       mjdn   = x_.nutdx.mjd;
       dX     = x_.nutdx.val;
       dX_e   = x_.nutdx.mx;
    else
       mjdnutdx   = [mjde(1);mjde(end)];
       dX     = [0; 0];
       dX_e   = [0; 0];
    end

    % nutation dY if estimated
    if opt_.nutdy.model==1
        mjdnutdy=x_.nutdy.mjd;
        dY     = x_.nutdy.val;
        dY_e   = x_.nutdy.mx;
    else
        mjdnutdy=[mjdxpol(1); mjdxpol(end)];
        dY     = [0; 0];
        dY_e   = [0; 0];
    end


    % interpolate all 
    xp   =interp1(mjdxpol,xp,mjd);
    xp_e =interp1(mjdxpol,xp_e,mjd);
    yp   =interp1(mjdypol,yp,mjd);
    yp_e =interp1(mjdypol,yp_e,mjd);
    ut   =interp1(mjddut1,ut,mjd);
    ut_e =interp1(mjddut1,ut_e,mjd);
    dX   =interp1(mjdnutdx,dX,mjd);
    dX_e =interp1(mjdnutdx,dX_e,mjd);
    dY   =interp1(mjdnutdy,dY,mjd);
    dY_e =interp1(mjdnutdy,dY_e,mjd);

elseif isempty(x_) && isempty(opt_) % if both are [] := normal without estimation   
    % set all esimates and errors to 0
    xp=zeros(size(mjd));
    yp=xp;
    ut=xp;
    dX=xp;
    dY=xp;
    xp_e=xp;
    yp_e=xp;
    ut_e=xp;
    dX_e=xp;
    dY_e=xp;
else
    % then both should exist (if only one is [] -> error message)
    disp('ERROR: Either x_ or opt_ does not exist!\n When nothing is estimated, both should be [], else both should exist!')
    keyboard;
end
% high frequency corrections
LEAP  = tai_utc(mjd);                   % get difference TAI-UTC [s] [25x1]
TT = mjd + (32.184 + LEAP')./86400;     % [MJD TT]  [1x25]
[xhf,yhf,uthf] = eophf(TT',parameter);   % [rad, sec]
xhf = rad2mas(xhf); yhf = rad2mas(yhf); uthf = uthf *1e3; %[mas, ms]  [25x1]each

% get a priori values
MJDeop =  parameter.eop.mjd;   %[12x1]
XPeop  = (parameter.eop.xp)*1e3;  % [mas]  [12x1]
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
    DUT1  = (lagint4v(MJDeop,UT1eop,mjd))';  %[1x25]
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
    XP   = XP;% + xhf';  %[1x25]
    YP   = YP;% + yhf';
    DUT1 = DUT1;% + uthf';
    inclhf = 'yes';
else
    inclhf = 'no';
end

if parameter.vie_mod.dXdY == 1
   inclanut = 'yes';
else
   inclanut = 'no';
end


% calculate total EOP values (tot = estimated + a priori + high frequency)
xptot = xp + XP;   % [mas]
yptot = yp + YP;   % [mas]
uttot = ut + DUT1; % [ms]
dXtot = dX + DXap; % [mas]
dYtot = dY + DYap; % [mas]


% create file
% unit: mas/ms
% out =
out = [mjd;xptot;yptot;uttot;dXtot;dYtot; XP;YP;DUT1;DXap;DYap; xp;yp;ut;dX;dY;xp_e;yp_e;ut_e;dX_e;dY_e;xhf';yhf';uthf'];

%      MJD|       TOTAL Values          |       apriori                 |    estimates |         errors         | high frequ 

%fid = fopen(strcat('../OUT/eop_',num2str(sname),'.txt'),'wt');
%fprintf('%%************************************************************\n');
%fprintf(strcat('%% EOP out table of lsm file x_',num2str(sname),'\n'));
% fprintf('%% Columns:\n%%\t 1     .... mjd\n%%\t 2-6   .... total values(x,y,ut,dX,dY)\n%%\t 7-11  .... a priori EOP (input in vie_mod)\n%%\t 12-16 .... estimated values\n%%\t 17-21 .... error of estimation\n%%\t 22-24 .... high frequency (subdaily) ERP corrections\n%%\n'); 
% fprintf('%% all units in mas resp. ms (dut1)\n');
% fprintf(strcat('%% subdaily high frequency variations included in the estimations: ',inclhf,'\n')); 
% fprintf(strcat('%% a priori nutation corrections were applied:',inclanut,'\n'));
% fprintf('%%************************************************************\n');
% fprintf('%% MJD\t\t xpol\t\t ypol\t\t dut1\t\t dX\t\t\t dY\t\t\t x_apr\t\t y_apr\t\t ut_apr\t\t dX_apr\t\t\t dY_apr\t\t x_est\t\t y_est\t dut1_est\t dX_est\t\t dY_est\t\t\t x_err\t\t y_err\t\t ut_err\t\t dX_err\t\t dY_err\t\t x_hf\t\t y_hf\t\t ut_hf\n');
% fprintf('%%\n');

%fprintf(' %5.4f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f % 11.6f\n',out);
%fclose(fid);

eop_struct=out;

%cd(curDir);

end