% *************************************************************************
%   Description:
%   Estimate Ionospheric Parameters
%
%   References: Hobiger, 2006: VLBI as a tool to probe the ionosphere
%               Dettmering et al., 2011: Systematic differences between
%               VTEC obtained by different space-geodetic techniques during
%               CONT08
%
%   Coded for VieVS: 
%   Oct 29 2012 by Benedikt Soja
%   
%  Revision 1: corrected small errors for better convergence in the
%  adjustment
%  Nov 05 2012 by Benedikt Soja
%  Revision 2: estimating two gradients now, Gn and Gs (see DGFI paper)
%  Nov 07 2012 by Benedikt Soja
%  Revision 3: corrected small error related to aposteriori sigmas
%  Mar 26 2013 by Benedikt Soja
%  - 2017-11-20, A. Hellerschmied: Call of function xyz2ell.m changed
% *************************************************************************

% *************************************************************************
%   Format of resulting structure "ion":
%
%   st_name: name of the station
%   mjd: time of the interval borders [modified julian date]
%   vtec: estimated vtec at the interval borders [TECU]
%   instr_bias: constant offset caused by instrumental delays [ns]
%   G_n/G_s: North/South gradient, usage: (1 + dlat*G_n) [1/rad]
%   *_sd: 1-sigma standard deviation of the estimated parateters
%
% *************************************************************************

% *************************************************************************
%   User Input
% *************************************************************************

%directory of your VieVS installation:
vievsdir = '../..';
%subfolder of the already processed session (can be an empty string):
subdir = '';
%name of the already processed session:
session = '11APR04XA_N004';

%scan and antenna mat-structures in DATA/LEVEL1 are necessary

%number of observations per interval
%lower values offer better temporal resolution, but may cause instabilities
%in the estimation
%values >= 8 (intervals of roughly 30min) are recommended
n_obs_int = 15; 

%maximum zenith distance in degrees
%all observations with higher zd are excluded from the estimation
%values <= 80deg are recommended
max_zd = 60; 

% *************************************************************************
%   End of User Input
% *************************************************************************


%% Loading and rearranging the data
clc;close all
clearvars -except vievsdir subdir session n_obs_int max_zd

%paths to mat-files 
path_s = strcat(vievsdir,'/DATA/LEVEL1/',subdir,'/',session,'_scan.mat');
path_a = strcat(vievsdir,'/DATA/LEVEL1/',subdir,'/',session,'_antenna.mat');

%load data or give warning
if exist(path_s,'file') ~= 2 || exist(path_a,'file') ~= 2 
    fprintf('\nSession is not yet processed or path name is wrong:\n\n')
    fprintf('%s\n',path_s)
    fprintf('%s\n\n',path_a)
    return
else
    load(path_s);
    load(path_a);
    fprintf('\nStarting processing of session %s...\n\n',session)
end


%constants
vx = 8.59549*10^9;  % X-band [Hz] 
% = mean([8.21299e+9 8.25299e+9 8.35299e+9 8.51299e+9 8.73299e+9 8.85299e+9 8.91299e+9 8.93299e+9])
c = 299792458; %speed of light [m/s]
factor = 40.28/c/vx^2*1e16*1e9; %factor 10^16 from TECU, 10^9 from delay in ns

%find out how many observations in total to initialize arrays
n=0;
for i = 1:length(scan) 
    n = n + length(scan(i).obs);
end

mjd=zeros(n,1); %time as modified julian date
i1=zeros(n,1); %station number
i2=zeros(n,1);
stc1=zeros(n,3); %station coordinates (x,y,z) in m 
stc2=zeros(n,3);
az1=zeros(n,1); %azimuth in rad
az2=zeros(n,1); 
zd1=zeros(n,1); %zenith distance in rad
zd2=zeros(n,1); 
dion=zeros(n,1); %ionospheric correction in ns
sdion=zeros(n,1); %std. dev. of ionospheric correction in ns
qdion=zeros(n,1); %quality code of ionospheric correction

%populate arrays with level 1 data
k=0;
for i = 1:length(scan)
    for j = 1:length(scan(i).obs)
        k=k+1;
        mjd(k) = scan(i).mjd;
        i1(k) = scan(i).obs(j).i1;
        i2(k) = scan(i).obs(j).i2;
        stc1(k,:) = scan(i).stat(i1(k)).x;
        stc2(k,:) = scan(i).stat(i2(k)).x;
        az1(k) = scan(i).stat(i1(k)).az;
        az2(k) = scan(i).stat(i2(k)).az;
        zd1(k) = scan(i).stat(i1(k)).zd;
        zd2(k) = scan(i).stat(i2(k)).zd;
        dion(k) = scan(i).obs(j).delion;
        sdion(k) = scan(i).obs(j).sgdion;
        qdion(k) = scan(i).obs(j).q_code_ion;
    end
end
clear scan

t = (mjd - floor(mjd(1)))*24; %time since 0 UT, day of session, in hours

n_st = max([i1;i2]); %number of stations

if sum(qdion) ~=0
    fprintf('Warning: quality code of some ionospheric corrections not zero\n');
end

%% Calculations necessary for setting up the adjustment

R = 6371; %earth radius in km
H=506.7; %height of single layer ionosphere in km
alpha=.9782; %factor in mf
% other possibility for mf (not recommended):
% H = 450; 
% alpha = 1; 

%ionospheric pierce point
[lat,lon,h] = xyz2ell(stc1);
stc1_ell = [lat,lon,h];
dpsi = zd1 - asin((R+stc1_ell(:,3)*1e-3)/(R+H).*sin(zd1));
dlat = asin(sin(dpsi).*cos(az1));
dlon = atan(tan(dpsi).*sin(az1))./cos(stc1_ell(:,1));

ipp1_ell = stc1_ell + [dlat dlon zeros(n,1)];
ipp1 = ell2xyz(ipp1_ell);

[lat,lon,h] = xyz2ell(stc2);
stc2_ell = [lat,lon,h];
dpsi = zd2 - asin((R+stc1_ell(:,3)*1e-3)/(R+H).*sin(zd2));
dlat = asin(sin(dpsi).*cos(az2));
dlon = atan(tan(dpsi).*sin(az2))./cos(stc2_ell(:,1));
ipp2_ell = stc2_ell + [dlat dlon zeros(n,1)];
ipp2 = ell2xyz(ipp2_ell);

%remove all observations where the ipp is across the pole
%or where the zd is too large
rem_ind = find(ipp1_ell(:,1)>pi/2 | ipp1_ell(:,1)<-pi/2 |...
    ipp2_ell(:,1)>pi/2 | ipp2_ell(:,1)<-pi/2 | ...
    zd1>max_zd*pi/180 | zd2>max_zd*pi/180 | sdion==0);

zd1(rem_ind) = [];
zd2(rem_ind) = [];
t(rem_ind) = [];
i1(rem_ind) = [];
i2(rem_ind) = [];
stc1(rem_ind,:) = [];
stc2(rem_ind,:) = [];
ipp1(rem_ind,:) = [];
ipp2(rem_ind,:) = [];
dion(rem_ind) = [];
sdion(rem_ind) = [];
n = n - length(rem_ind);


%transformation of station and ipp coordinates to a geomagnetic system
%cartesian geomagnetic
stc1m = geogr2magn(stc1);
stc2m = geogr2magn(stc2);
ipp1m = geogr2magn(ipp1);
ipp2m = geogr2magn(ipp2);
%lat,lon geomagnetic
[lat,lon,h] = xyz2ell(stc1m);
stc1me = [lat,lon,h];
[lat,lon,h] = xyz2ell(stc2m);
stc2me = [lat,lon,h];
[lat,lon,h] = xyz2ell(ipp1m);
ipp1me = [lat,lon,h];
[lat,lon,h] = xyz2ell(ipp2m);
ipp2me = [lat,lon,h];


%"corrected" observation time  with longitude difference
dlon1 = ipp1me(:,2) - stc1me(:,2);
dlon2 = ipp2me(:,2) - stc2me(:,2);
%remove jumps of 2*pi, assumption that max_zd is <= 85°, Ny Alesund at 79°N
dlon1 = dlon1 + (dlon1<-pi/2)*2*pi - (dlon1>pi/2)*2*pi;
dlon2 = dlon2 + (dlon2<-pi/2)*2*pi - (dlon2>pi/2)*2*pi;
t1 = t + dlon1*180/pi/15; %in h
t2 = t + dlon2*180/pi/15; 

%latitude difference
dlat1 = ipp1me(:,1) - stc1me(:,1); %in rad
dlat2 = ipp2me(:,1) - stc2me(:,1);

%indeces for observations with N or S gradients
igr1 = [dlat1 >= 0 dlat1 < 0];
igr2 = [dlat2 >= 0 dlat2 < 0];

%mapping function
m1 = 1./sqrt(1 - (R/(R+H)*sin(alpha*zd1)).^2);
m2 = 1./sqrt(1 - (R/(R+H)*sin(alpha*zd2)).^2);


%% Design matrix and adjustment

%find number of observations and intervals for each station
imat = -(i1(:,ones(1,n_st))==ones(n,1)*(1:n_st)) + ...
    (i2(:,ones(1,n_st))==ones(n,1)*(1:n_st)); %index-matrix (1 for stat2, -1 for stat1)
n_stobs = sum(abs(imat)); %number of observations per station

n_int = floor(n_stobs/n_obs_int); %number of intervals per station
i_unkn = [1 1+4*(1:n_st-1)+cumsum(n_int(1:end-1))]; %index of first unknown per station

%catch error for not enough obs:
if sum(n_int == 0)
    fprintf('Not enough observations available for the following stations:\n\n')
    for i = find(n_int == 0)
        fprintf('%s\n',antenna(i).name);
    end
    fprintf('\nProcessing cancelled!\n\n')
    return
end

%start and end times for each interval
t_all = cell(1,n_st);
t_int = nan(max(n_int)+1,n_st);
for i=1:n_st %loop over stations
    t_all{i} = sort([t1(imat(:,i)==-ones(n,1)); t2(imat(:,i)==ones(n,1))]);
    t_int(1:n_int(i)+1,i) = [t_all{i}(1); mean([t_all{i}((1:n_int(i)-1)*n_obs_int)...
        t_all{i}((1:n_int(i)-1)*n_obs_int+1)],2); t_all{i}(end)];
end

%index of the last interval for each observation
i_rt1 = sum(bsxfun(@lt,t_int(:,i1)',t1),2);
i_rt2 = sum(bsxfun(@lt,t_int(:,i2)',t2),2);

%differences between intervals for all observations
int_diffs1 = t_int(2:end,i1)' - t_int(1:end-1,i1)';
int_diffs2 = t_int(2:end,i2)' - t_int(1:end-1,i2)';
for i=1:n
    if i_rt1(i) == 0
        int_diffs1(i,:) = zeros(1,max(n_int));
    else        
        int_diffs1(i,i_rt1(i):end) = [t1(i)-t_int(i_rt1(i),i1(i)) zeros(1,max(n_int)-i_rt1(i))];
    end
    if i_rt2(i) == 0
        int_diffs2(i,:) = zeros(1,max(n_int));
    else        
        int_diffs2(i,i_rt2(i):end) = [t2(i)-t_int(i_rt2(i),i2(i)) zeros(1,max(n_int)-i_rt2(i))];
    end
end


%B matrix (non-negative constraints)
B = zeros(n_st+sum(n_int),4*n_st+sum(n_int));
for i = 1:n_st
   B(i_unkn(i)-(i-1)*2:i_unkn(i)-(i-1)*2+n_int(i), i_unkn(i)+3) = ones(n_int(i)+1,1);
   B(i_unkn(i)-(i-1)*2+1:i_unkn(i)-(i-1)*2+n_int(i), i_unkn(i)+4:i_unkn(i)+4+n_int(i)-1) = ...
       triu(toeplitz(ones(1,n_int(i))))'*diag((t_int(2:n_int(i)+1,i)' - t_int(1:n_int(i),i)'));
end

%P-Matrix
k=6;
wk = (2*(1-(R/(R+H))^2*sin(zd2).^2).*(1-(R/(R+H))^2*sin(zd1).^2)./...
    (2-(R/(R+H))^2*(sin(zd2).^2+sin(zd1).^2))).^(k/2);
P = diag([wk./sdion.^2; 1]); %1 for datum row

A = zeros(n+1,4*n_st+sum(n_int));
A(end,i_unkn)=ones(1,n_st); %datum of instrumental biases: their sum is zero
x0 = zeros(4*n_st+sum(n_int),1); %initial guess for unknowns, used for L0, A
x0(i_unkn+1) = 1; %Gn
x0(i_unkn+2) = 1; %Gs
x0(i_unkn+3) = 20; %offs
L0 = zeros(n+1,1); %computed ionospheric delays
s0_ = 100; %s0 of the "previous iteration"

dt_times_rt1 = zeros(n,1);
dt_times_rt2 = zeros(n,1);

%options for quadprog
options = optimset('Display','off','MaxIter',20);
warning off all

max_iter = 15;
for k=1:max_iter %iterative adjustment
    
    %L0
    for i = 1:n
        dt_times_rt1(i)=int_diffs1(i,1:n_int(i1(i)))*x0(i_unkn(i1(i))+4:i_unkn(i1(i))+n_int(i1(i))+3);
        dt_times_rt2(i)=int_diffs2(i,1:n_int(i2(i)))*x0(i_unkn(i2(i))+4:i_unkn(i2(i))+n_int(i2(i))+3);
    end
    L0 = factor*(-m1.*(1 + sum(igr1.*[x0(i_unkn(i1)+1) x0(i_unkn(i1)+2)],2).*dlat1).*...
        (x0(i_unkn(i1)+3) + dt_times_rt1) +...
        m2.*(1 + sum(igr2.*[x0(i_unkn(i2)+1) x0(i_unkn(i2)+2)],2).*dlat2).*...
        (x0(i_unkn(i2)+3) + dt_times_rt2)) - x0(i_unkn(i1)) + x0(i_unkn(i2));
    
    %obs-comp, 0 for datum (sum of instrumental biases = 0)
    l = [dion - L0; 0];
    
    %Design matrix A
    for i = 1:n
       A(i,i_unkn(i1(i)):i_unkn(i1(i))+n_int(i1(i))+3) = -1*[1 ... %instr bias
           factor*m1(i)*igr1(i,:)*dlat1(i)*(x0(i_unkn(i1(i))+3) + dt_times_rt1(i))... %gradients
           factor*m1(i)*(1+(igr1(i,:)*[x0(i_unkn(i1(i))+1);x0(i_unkn(i1(i))+2)])*dlat1(i))... offset
           factor*m1(i)*(1+(igr1(i,:)*[x0(i_unkn(i1(i))+1);x0(i_unkn(i1(i))+2)])*dlat1(i))*int_diffs1(i,1:n_int(i1(i)))]; %rates
       
       A(i,i_unkn(i2(i)):i_unkn(i2(i))+n_int(i2(i))+3) = [1 ...
           factor*m2(i)*igr2(i,:)*dlat2(i)*(x0(i_unkn(i2(i))+3) + dt_times_rt2(i))...
           factor*m2(i)*(1+(igr2(i,:)*[x0(i_unkn(i2(i))+1);x0(i_unkn(i2(i))+2)])*dlat2(i))...
           factor*m2(i)*(1+(igr2(i,:)*[x0(i_unkn(i2(i))+1);x0(i_unkn(i2(i))+2)])*dlat2(i))*int_diffs2(i,1:n_int(i2(i)))];
    end
    x_ = quadprog(A'*P*A,-(A'*P*l)',-B,B*x0,[],[],[],[],[],options); %adjustment
    x0 = x0 + x_;
    v = A*x_- l; 
    s0 = sqrt(v'*P*v/(length(dion)+1-length(x_))); %std.dev. of unit weight
    fprintf('Iteration %d: s0 = %4.3g cm\n',k-1,s0*c*1e-9*100)

    %stop if change in s0 < 0.1mm
    if abs(s0 - s0_)*c*1e-9*100 < 0.01
        break
    else
        s0_ = s0;
    end
end

if k == max_iter
    fprintf('\nWarning: maximum number of iterations reached!\n');
end

fprintf('\nProcessing finished. Results can be found in the "ion" structure.\n')
    

 %% extracting variables and plotting results

Qxxd = diag(inv(A'*P*A));
t_min = min(min(t_int));
t_max = max(max(t_int));
tec_max = max(max(B*x0));

ion = struct('st_name',{},'mjd',{},'vtec',{},'vtec_sd',{},'instr_bias',{},...
    'instr_bias_sd',{},'G_n',{},'G_n_sd',{},'G_s',{},'G_s_sd',{});

figure
for i=1:n_st
    ind = i_unkn(i)-(i-1)*2:i_unkn(i)-(i-1)*2+n_int(i);
    subplot(ceil(n_st/2),2,i)
    plot(t_int(1:n_int(i)+1,i),B(ind,:)*x0,'.-')
    hold on
    title(antenna(i).name)
    ylabel('VTEC [TECU]')
    if i == n_st-1 || i == n_st
        xlabel('t since 0 UT, day 1 [h]')
    end
    axis([t_min t_max 0 tec_max])
    
    var = Qxxd(i_unkn(i):i_unkn(i)+3+n_int(i))*s0^2; %variance of unknowns
    %populate ion struct
    ion(i).st_name = antenna(i).name;
    ion(i).mjd = t_int(1:n_int(i)+1,i)/24 + floor(mjd(1));
    ion(i).vtec = B(ind,:)*x0;
    ion(i).instr_bias = x0(i_unkn(i));
    ion(i).G_n = x0(i_unkn(i)+1);
    ion(i).G_s = x0(i_unkn(i)+2);
    ion(i).instr_bias_sd = sqrt(var(1));
    ion(i).G_n_sd = sqrt(var(2));
    ion(i).G_s_sd = sqrt(var(3));
    %error propagation (from offs+rt to vtec):
    ion(i).vtec_sd = sqrt(B(ind,i_unkn(i):i_unkn(i)+3+n_int(i)).^2*var);
end
clearvars -except ion
warning on all
