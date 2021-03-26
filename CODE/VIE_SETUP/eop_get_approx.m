% #########################################################################
% #     eop_get_approx
% #########################################################################
%
% DESCRIPTION
%   EOP approximation of a priori
% 
%
% CREATED  
%   2016/11     A. Girdiuk
%
% INPUT
% - parameter   : structure from LEVEL3
% - mjd         : vector of time momemts of observations 
% - TT          : mjd vector converted to TT
%
% OUTPUT
%       (XP, YP) -polar motion, DUT1 -changes in UT1, 
%       (DXap, DYap) - precesion and nutation, 
%        MJDeop - corresponded time vector
%
% CHANGES:
%

function [XP,YP,DUT1,DXap,DYap,MJDeop] = eop_get_approx(parameter,mjd,TT)

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
%    disp('remove tidal UT')
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
%    disp('re-add tidal UT')
    %corrUT1 = tver2000(TT);     % [sec]
    corrUT1 = rg_zont2(TT,par35);  % [sec]
    DUT1    = DUT1 + corrUT1'*1e3;   % [ms]
end

%	if parameter.vie_mod.eophf ==1;
%		inclhf = 'yes';
%	else
%		inclhf = 'no';
%	end

%	if parameter.vie_mod.dXdY == 1
%	   inclanut = 'yes';
%	else
%	   inclanut = 'no';
%	end


