% ************************************************************************
%   Description:
%   This function first generates a vector of turbulent equivalent zenith
%   wet delays (ezwd) for a certain station. The ezwd is then mapped to the
%   elevation of the observation using VMF1 or GMF. The mapping function is
%   supposed to be free of error (i.e. the same mapping function and
%   coefficients are used for the generation of the slant delays as are for
%   the lsm adjustment).
%
%   References: 
%    - Nilsson et al. 2007 "Simulations of atmospheric path delays using
%       turbulence models"
%    - Böhm et al. 2007 "Simulation of zenith wet delays and clocks"
%    - Treuhaft and Lanyi 1987 "The effect of the dynamic wet troposphere
%      on radio interferometric measurements"
%
%   Input:		
%      mjd:        [d] vector time in mjd
%      az:         [rad] vector azimuth
%      el:         [rad] vector elevation
%      mfw:        [ ] mapping function wet
%      idays:      [d] how many days should be simulated
%      dhseg:      [h] time segments over which observations are to be correlated in
%                  hours. if 2 hours is not short enough it might be better to
%                  use 1. If possible use larger values. 24 hours would be
%      dh:         [m] height increment for numeric integration (e.g. 200 m)
%                  optimum, but perhaps not feasible.
%      Cn:         [1e-7*m^(-1/3)] refractive index structure constant
%      h:          [m] effective height of the troposphere
%      vn:         [m/s] north component of the wind vector
%      ve:         [m/s] east vomponent of the wind vector
%      wzd0:       [mm] initial zenith wet delay
% 
%   Output:
%      a vector of simulated equivalent zenith wet delays [mm]
% 
%   External calls: 	
%      constants.m   
%       
%   Coded for VieVS: 
%   adopted to VieVS July 2010 by Andrea Pany
%   (original function coded by J. Böhm, Sep 6 2007)
%
%   Revision: 
%   ---
% ************************************************************************

function [swd,cov] = sim_swd(mjd,az,el,mfw,idays,dhseg,dh,Cn,h,vn,ve,wzd0)

% constants and settings
L0    = 3e6;               % [m] saturation scale length (Treuhaft and Lanyi 1987)
Cn    = Cn*1e-7;           % [1e-7*m^(-1/3)]              
zwd   = [];

mjd = mjd(:); az = az(:); el = el(:);

% constants: light velocity     % Jing SUN, Jan 11, 2012
c = 299792458;   %[m/s]         % Jing SUN, Jan 11, 2012

% open the masterfile for the stations to get the information
vn = vn*3600;           % [m/h] north component of wind vector
ve = ve*3600;           % [m/h] east component of wind vector
vv = [vn ve 0];         % wind vector

% number of observations for the specific station
num = length(mjd);

% time epochs of the observations in hours
t = mjd - floor(mjd);
for i = 2:num
    if t(i) < t(i-1)
        t(i) = t(i) + 1;
    end
end
t = t*24;

ri = zeros(num,3);
for i = 1:num
    ri(i,1:3) = [cos(az(i))/tan(el(i)) sin(az(i))/tan(el(i)) 1];
end
% ---
L23 = L0^(2/3);
Cnall = Cn^2/2*1e6*dh^2;
% number of height increments
nh    = floor(h/dh);       
% the first observation is sacrificed :-)
az1 = 0;
el1 = pi/2;
r0 = [cos(az1)/tan(el1) sin(az1)/tan(el1) 1];
k = 0;
% preallocation
zs = zeros(1,(nh+1)^2);
z  = zeros(1,(nh+1)^2);
v  = zeros(3,(nh+1)^2);
for i = 1:nh+1
    for j = 1:nh+1
        k = k + 1;
        zs(k) = dh*(i-1);
        z (k) = dh*(j-1);
        v(1:3,k) = vv';
    end
end
r0z  = r0'*z;
r0zs = r0'*zs;
rho4  = sqrt(sum((r0z-r0zs).^2));
rho4x = rho4.^(2/3)./(1+rho4.^(2/3)/L0^(2/3));
% --
% split everything into two-hour segments ...
ti = floor(t./dhseg)+1;

% how many epochs per segment
isegments = ti(num);
% preallocation
tn = zeros(1,isegments);
for i = 1:isegments
    k = find(ti==i);
    tn(i) = length(k);
end
% if only one segment
if isegments==1
    num3 = tn(1); 
    C(1:num3,1:num3) = 0;
    % how many epochs are before the segment i1
    for i = 1:num3
        riz  = ri(i,1:3)'*z;
        dti0 = t(i);
        dd1 = (riz-r0zs+v*dti0);
        rho1  = sum(dd1.*dd1).^(1/3);
        rho1x = rho1./(1+rho1/L23);
        for j = i:num3
            % compute time differences in hours
            dtj0 = t(j);
            dtij = t(j)-t(i);
            % compute position vectors
            rjz  = ri(j,1:3)'*z;
            rjzs = ri(j,1:3)'*zs;
            % compute separations
            dd2 = (rjz-r0zs+v*dtj0);
            dd3 = (riz-rjzs-v*dtij);
            rho2 = sum(dd2.*dd2).^(1/3);
            rho3 = sum(dd3.*dd3).^(1/3);
            out  = rho1x + rho2./(1+rho2/L23) - rho3./(1+rho3/L23) - rho4x;%
            result = sum(out);
            C(i,j) = Cnall*result;
        end
    end
    D = chol(C)';
    cov = D;
    swd = zeros(num,idays);
    for iday = 1:idays
        x = randn(num3,1);
        l = D*x;
        % add an initial value
        swd(:,iday) =((l + wzd0)'.*mfw)';
    end
else
% --
    % first two segments
    i1 = 1;
    num1 = tn(i1);
    num2 = tn(i1+1);
    num3 = num1+num2;
%     C = []; D = []; D11 = []; D12 = []; D22 = [];
    C = zeros(num3,num3);
    % how many epochs are before the segment i1
    k = sum(tn(1:i1-1));
    for i = 1:num3
        riz  = ri(k+i,1:3)'*z;
        dti0 = t(k+i);
        dd1 = (riz-r0zs+v*dti0);
        rho1  = sum(dd1.*dd1).^(1/3);
        rho1x = rho1./(1+rho1/L23);
        for j = i:num3
            % compute time differences in hours
            dtj0 = t(k+j);
            dtij = t(k+j)-t(k+i);
            % compute position vectors
            rjz  = ri(k+j,1:3)'*z;
            rjzs = ri(k+j,1:3)'*zs;
            % compute separations
            dd2 = (rjz-r0zs+v*dtj0);
            dd3 = (riz-rjzs-v*dtij);
            rho2 = sum(dd2.*dd2).^(1/3);
            rho3 = sum(dd3.*dd3).^(1/3);
            out  = rho1x + rho2./(1+rho2/L23) - rho3./(1+rho3/L23) - rho4x;%
            result = sum(out);
            C(i,j) = Cnall*result;
        end
    end
    D = chol(C)';
    C11 = C(num1+1:num3,num1+1:num3);
    cov(i1).D = sparse(D);
    cov(i1).num1 = num1;
    cov(i1).num2 = num2;
    cov(i1).num3 = num3;
    % --
    % loop over two adjacent two-hour blocks
    for i1 = 2:isegments-1
        num1 = tn(i1);
        num2 = tn(i1+1);
        num3 = num1+num2;
%         C = []; D = []; D11 = []; D12 = []; D22 = [];
        C = zeros(num3,num3);
        % how many epochs are before the segment i1
        k = sum(tn(1:i1-1));
        C(1:num1,1:num1) = C11;
        for i = 1:num1
            riz  = ri(k+i,1:3)'*z;
            dti0 = t(k+i);
            dd1 = (riz-r0zs+v*dti0);
            rho1  = sum(dd1.*dd1).^(1/3);
            rho1x = rho1./(1+rho1/L23);
            for j = num1+1:num3
                % compute time differences in hours
                dtj0 = t(k+j);
                dtij = t(k+j)-t(k+i);
                % compute position vectors
                rjz  = ri(k+j,1:3)'*z;
                rjzs = ri(k+j,1:3)'*zs;
                % compute separations
                dd2 = (rjz-r0zs+v*dtj0);
                dd3 = (riz-rjzs-v*dtij);
                rho2 = sum(dd2.*dd2).^(1/3);
                rho3 = sum(dd3.*dd3).^(1/3);
                out  = rho1x + rho2./(1+rho2/L23) - rho3./(1+rho3/L23) - rho4x;
                result = sum(out);
                C(i,j) = Cnall*result;
            end
        end
        for i = num1+1:num3
            riz  = ri(k+i,1:3)'*z;
            dti0 = t(k+i);
            dd1 = (riz-r0zs+v*dti0);
            rho1  = sum(dd1.*dd1).^(1/3);
            rho1x = rho1./(1+rho1/L23);
            for j = i:num3
                % compute time differences in hours
                dtj0 = t(k+j);
                dtij = t(k+j)-t(k+i);
                % compute position vectors
                rjz  = ri(k+j,1:3)'*z;
                rjzs = ri(k+j,1:3)'*zs;
                % compute separations
                dd2 = (rjz-r0zs+v*dtj0);
                dd3 = (riz-rjzs-v*dtij);
                rho2 = sum(dd2.*dd2).^(1/3);
                rho3 = sum(dd3.*dd3).^(1/3);
                out  = rho1x + rho2./(1+rho2/L23) - rho3./(1+rho3/L23) - rho4x;
                result = sum(out);
                C(i,j) = Cnall*result;
            end
        end
        D = chol(C)';
        C11 = C(num1+1:num3,num1+1:num3);
        cov(i1).D = sparse(D);
        cov(i1).num1 = num1;
        cov(i1).num2 = num2;
        cov(i1).num3 = num3;
    end
    % how many days should be simulated
    swd = zeros(num,idays);
    for iday = 1:idays
        x = randn(cov(1).num3,1);
        l = cov(1).D*x;
        l1 = l(cov(1).num1+1:cov(1).num3);
        for i1 = 2:isegments-1
            D11 = cov(i1).D(1:cov(i1).num1,1:cov(i1).num1);
            D21 = cov(i1).D(cov(i1).num1+1:cov(i1).num3,1:cov(i1).num1);
            D22 = cov(i1).D(cov(i1).num1+1:cov(i1).num3,cov(i1).num1+1:cov(i1).num3);
            x = randn(cov(i1).num2,1);
            l1 = D21*inv(D11)*l1 + D22*x;
            l = [l; l1];
        end
        % add an initial value
        swd(:,iday) = ((l + wzd0)'.*mfw)';
    end
end

% % map equivalent zenith wet delay to the elevations of the observations
% swd = zwd'.*mfw;      % [mm]
% convert swd from [mm] to [s]
swd = swd*1e-3*1/c;