% ************************************************************************
%   Description:
%   function to interpolate EOP to observation epochs
%
%   Reference:
%
%   Input:
%       'parameter'      structure array    parameter structure
%       'MJDeop'         (n,1)              MJD loaded from a-priori time series
%       'UT1eop'         (n,1)              UT1 loaded from a-priori time series
%       'XPeop'          (n,1)              XP loaded from a-priori time series
%       'YPeop'          (n,1)              YP loaded from a-priori time series
%       'dXeop'          (n,1)              dX loaded from a-priori time series
%       'dYeop'          (n,1)              dY loaded from a-priori time series
%       'MJD'            (x,1)              MJD time of observations       
%       'TT'             (x,1)              Terrestrial time of observations
%       'eopch_type'     string             type of obervations - either observation (quasar) or satellite (orbit)
%
%   Output:
%       'parameter'      structure array    sources structure          
%       'DUT1'           (x,1)              DUT1 at obervation epochs
%       'XP'             (x,1)              XP at observation epochs
%       'YP'             (x,1)              YP at obsevation epochs
%       'DX'             (x,1)              DX at observation epochs
%       'DY'             (x,1)              DY at observation epochs
%
%   External calls:
%       rg_zont2.m   lagint4v.m   eophf.m
% 
%
%   Coded for VieVS:
%   22 November 2021 by H. Wolf - created as external function of vie_mod
%
%   Revision:
%
% ************************************************************************


function [parameter,DUT1, XP, YP, DX, DY] = interpolateEOP(parameter, MJDeop, UT1eop, XPeop, YPeop, dXeop, dYeop, MJD, TT, epoch_type)
    % subtraction of tidal variations (Defraigne and Smits) in dUT1 before interpolation:
    if parameter.vie_mod.tidalUT == 1
        taiut    = tai_utc(MJDeop);
        MJDTTeop = MJDeop + (32.184 + taiut)/86400;
        %         UT1corr  = tver2000(MJDTTeop);  % [sec]
        if parameter.vie_mod.tidalUT35 ==1
            par35=1;
        else
            par35=2;
        end
        UT1corr  = rg_zont2(MJDTTeop, par35);     % [sec]
        UT1eop   = UT1eop - UT1corr;
    end

    % Interpolation:
    % Linear
    if parameter.vie_mod.linear == 1
        disp('linear interpolation of EOP')
        parameter.eop.interp = 'linear';
        if parameter.vie_mod.linear48h && strcmp(epoch_type,'observation')
            disp('special linear interpolation for IVS SINEX submission with 48h EOP interval')
            % find index of MJDeop that lies within the 48 hour estimation interval
            midof48h = find(MJDeop==floor(min(MJD))+1);
            % remove this midnight point from a priori time series
            MJDeop(midof48h) = []; UT1eop(midof48h) = [];
            XPeop(midof48h)  = []; YPeop(midof48h)  = []; 
            dXeop(midof48h)  = []; dYeop(midof48h)  = [];

            % store a priori values at EOP timesteps without the 24h midnight
            parameter.eop.mjd(midof48h) =  [];
            parameter.eop.xp(midof48h)  =  [];
            parameter.eop.yp(midof48h)  =  [];
            parameter.eop.ut1(midof48h) =  [];
            parameter.eop.dX(midof48h)  =  [];
            parameter.eop.dY(midof48h)  =  [];
        end
        % A priori EOP values are determined with linear interpolation between
        % the value of midnight before and after observation time.
        % for a session from 18:00 to 18:00 this means, that there are 2 a
        % priori lines and a break at midnight
        % no a priori values are stored in the parameter file!!!
        DUT1  = interp1(MJDeop,UT1eop,MJD,'linear','extrap');
        XP    = interp1(MJDeop, XPeop,MJD,'linear','extrap');
        YP    = interp1(MJDeop, YPeop,MJD,'linear','extrap');
        DX    = interp1(MJDeop, dXeop,MJD,'linear','extrap');
        DY    = interp1(MJDeop, dYeop,MJD,'linear','extrap');

    % Lagragne
    else % linear = 0
        disp('Lagrange interpolation of EOP')
        parameter.eop.interp = 'lagrange';
        % subtraction of tidal variations (Defraigne and Smits) in dUT1 before interpolation
        % interpolate EOP for time of observation
        DUT1    = lagint4v(MJDeop, UT1eop,MJD);
        XP      = lagint4v(MJDeop, XPeop, MJD);
        YP      = lagint4v(MJDeop, YPeop, MJD);
        DX      = lagint4v(MJDeop, dXeop, MJD);
        DY      = lagint4v(MJDeop, dYeop, MJD);
    end

    % re-add tidal variation in dUT1 after interpolation:
    if parameter.vie_mod.tidalUT == 1
        corrUT1 = rg_zont2(TT,par35);  % [sec]
        DUT1    = DUT1 + corrUT1;      % [sec]
    end

    % Add high frequency EOP
    [CORX,CORY,CORUT,parameter] = eophf(TT, parameter);
    DUT1 = DUT1 + CORUT;
    XP   = XP + CORX;
    YP   = YP + CORY;
end

