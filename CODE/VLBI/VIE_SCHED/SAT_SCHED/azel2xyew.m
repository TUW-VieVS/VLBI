% #########################################################################
% #     azel2xyew
% #########################################################################
%
% DESCRIPTION
%   Conversion of Az/El antenna pointing angles to XY- angels of XYEW antenna mount. 
%   Also change rates are calculated.
%
% CREATED  
%   2015-06-23     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%
%
% INPUT
% - az              - Azimuth, [deg] or [rad]
% - el              - Elevartion, , [deg] or [rad]
% - flag_in_degrees - If this flag is set (=1), output and input values are assumed to be in [degrees]! Otherwise in [rad].
% - az_rate         - Azimuth rate, [deg/sec] or [rad/sec]
% - el_rate         - Elevartion, [deg/sec] or [rad/sec]
%
%
% OUTPUT
% - x               - X angle [deg] or [rad]
% - y               - Y angle [deg] or [rad]
% - error_code      - Error Code (0 = no erros occured)
% - error_msg       - Error Message (empty, if no errors occured)
% - x_rate          - X angle rate [deg/sec] or [rad/sec]
% - y_rate          - Y angle rate [deg/sec] or [rad/sec]
%
% CHANGES:
%
%
function [x, y, error_code, error_msg, x_rate, y_rate] = azel2xyew(az, el, varargin)

    % init:
    error_code = 0;
    error_msg = '';
    az_rate = [];
    el_rate = [];
    x = [];
    y = [];
    x_rate = [];
    y_rate = [];
    flag_calc_rates = 0;
    flag_in_degrees = 0; % default = [rad]

    % ##### Check input #####
    switch(nargin)
        
        case 2 % In: only az/el;                            out: x/y [rad]
            
        case 3 % In: az/el, flag_in_degrees;                out: x/y in [rad] or [deg]
            flag_in_degrees = varargin{1};
            
        case 5 % In: az/el, flag_in_degrees, az/el-rates;   out: x/y, x/y-rates in [rad] or [deg]
            flag_in_degrees = varargin{1};
            az_rate         = varargin{2};
            el_rate         = varargin{3};
            flag_calc_rates = 1;
            
        otherwise % Error
            error_code = 1;
            error_msg = 'Incorrect input!';
            return;
            
    end % switch(nargin)
    
    
    % ##### Conversion: Az/El => X/Y_ew #####
    % formulars provided by Jamie McCallum (June 2015).
    
    if flag_in_degrees
        % Conversion [deg] => [rad]:
        az = az *pi/180;
        el = el *pi/180;
    end
    
        
    % Code copied from:
    % --------------------------------------------------------------------
    %   function [xcor,ycor, x30, y30] = aziel2xy(azi,ele)
    %   % Converts Azimuth and elevation to X Y coordinates
    %   % Written by Cliff Senkbeil (2004)
    n = length(az);
    m = length(el);
    if n ~= m
        error('Azimuth dimension not equal to Elevation dimension')
    end
    polax = cos(az).*cos(el);
    polay = sin(az).*cos(el);
    polaz = sin(el);
    for jpola = 1:n
        theta(jpola)=atan2(polax(jpola),polaz(jpola));
        radius(jpola)=sqrt(polax(jpola)^2+polay(jpola)^2+polaz(jpola)^2); % is always = 1!!!
        phi(jpola)=asin(polay(jpola)/radius(jpola));
%         theta(jpola);
%         phi(jpola);
    end
    for i = 1:n
        if theta(i) > pi
            theta(i) = theta(i) - 2*pi;
        end
        if phi(i) > pi
            phi(i) = phi(i) - 2*pi;
        end
    end
    x=theta';
    y=phi';
    % --------------------------------------------------------------------
    
    if flag_in_degrees
        % Conversion [rad] => [deg]:
        x = x /pi *180;
        y = y /pi *180;
    end
    
    
    
    % ##### Calculation of the x/y-rates #####
    if flag_calc_rates

        if flag_in_degrees
            % Conversion [deg] => [rad]:
            az_rate = az_rate *pi/180;
            el_rate = el_rate *pi/180;
        end
        
        saz = sin(az);
        sel = sin(el);
        cel = cos(el);
        caz = cos(az);
       
        % y_rate [rad/sec]:
        y_rate = 1/sqrt(1 - (saz*cel)) * (caz*cel*az_rate - saz*sel*el_rate);
        
        % x_rate [rad/sec]:
        x_rate = 1/(1 + ((caz*cel)/sel)^2) * -(((saz*cel*az_rate + caz*sel*el_rate)/sel) + ((caz*cel^2*el_rate)/sel^2)); 
        
        if flag_in_degrees
            % Conversion [rad/sec] => [deg/sec]:
            y_rate = y_rate /pi*180;
            x_rate = x_rate /pi*180;
        end
        
    end
    
return
