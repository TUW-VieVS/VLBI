% ************************************************************************
%   Description:
%   Solid Earth Tides corrections; displacement caused by lunar and solar
%   gravitational attraction.
% 
%   Reference: 
%   IERS Conventions 2010, Chapter 7.1.1
%
%   Input:										
%      mjd                 Modified Julian Date [d]
%      leap                difference between UTC and TAI in [s]
%                          leap seconds(TAI-UTC))
%      t2c                 transformation matrix TRS --> CRS [3,3]
%      ant                 station cartesian coordinates from catalogue, 
%                          TRF [m]
%      xyz_moon            CRS coordinates of the Moon at obs. time [m]
%      xyz_sun             CRS coordinates of the Sun  at obs. time [m]
%
%   Output:
%      cts                 station displacement vector (x,y,z) [m]
% 
%   External calls: 	
%      global rearthm massrSE massrME
%      doodarg.m, cart2phigd.m, ren2xyz.m
%
%   Coded for VieVS: 
%   19 Jun 2012 by Hana Krásná (function following the IERS 2010 Conventions (without any partial derivatives))
%

% *************************************************************************
%
function [cts] = mathews(mjd,leap,t2c,ant,xyz_moon,xyz_sun)
  
global rearthm massrSE massrME 

% factors
Re  = rearthm;      %[m]
MRs = massrSE;      %mass ratio Sun/Earth
MRm = massrME;      %mass ratio Moon/Earth


% Love numbers
h_0 = 0.6078;     %h^(0)
h_2 =-0.0006;      %h^(2)
h3  = 0.292;

l_0 = 0.0847;     %l^(0)
l_2 = 0.0002;     %l^(2)
l3  = 0.015;

l_11 = 0.0012;     %l^(1) diurnal band
l_12 = 0.0024;     %l^(1) semidiurnal band 

hI_21=-0.0025;
lI_21=-0.0007;
hI_22=-0.0022;
lI_22=-0.0007;

% unit vectors REN
r=[1,0,0];
e=[0,1,0];
n=[0,0,1];

%%

% rotate CRS coordinate to TRS 
moon = t2c'*xyz_moon;
sun  = t2c'*xyz_sun;

%geocentric east longitude, latitude, distance Moon/Sun
[Lm,Bm,Rm] = cart2sph(moon(1),moon(2),moon(3));          %[rad,rad,m]
[Ls,Bs,Rs] = cart2sph(sun(1),sun(2),sun(3));
%longitude,latitude,Ra for antenna
[lam,phi,Ra] = cart2sph(ant(1),ant(2),ant(3)); 

%unit vector from the geocentre to the Moon/Sun/antenna
rm = moon/Rm;
rs = sun/Rs;
ra = ant/Ra;

% scalar product of unit vectors antenna and Moon/Sun
scmon=dot(ra,rm);
scsun=dot(ra,rs);

% factors
f2m=MRm*Re^4/Rm^3;
f2s=MRs*Re^4/Rs^3;
f3m=f2m*Re/Rm;
f3s=f2s*Re/Rs;

%Legendre polynomial of degree 2, for sin(phi), latitude of antenna
P2_sphi=3/2*(sin(phi))^2-1/2;

% cm: cosinus of angle between antenna and Moon
% cs: cosinus of angle between antenna and Sun
P2_cm=3/2*scmon^2-1/2;
P2_cs=3/2*scsun^2-1/2;

%Legendre polynomial of degree 3
P3_cm=5/2*scmon^3-3/2*scmon;
P3_cs=5/2*scsun^3-3/2*scsun;

%Legendre function of degree 2, m=1,2
P21_sBm=3*cos(Bm)*sin(Bm);
P21_sBs=3*cos(Bs)*sin(Bs);
P22_sBm=3*(cos(Bm))^2;
P22_sBs=3*(cos(Bs))^2;

% new h2 and l2, with station latitude phi
h2=h_0+h_2*P2_sphi;
l2=l_0+l_2*P2_sphi;

%%
% IERS10 eq.7.5
% Displacement due to degree-2 tides, with nominal values for h2 and l2
dx9_m=f2m*(h2*ra*P2_cm+3*l2*scmon*(rm-scmon*ra')'); %m
dx9_s=f2s*(h2*ra*P2_cs+3*l2*scsun*(rs-scsun*ra')'); %m
dx9=dx9_m+dx9_s; %m

% IERS10 eq.7.6
% Displacement due to degree-3 tides
dx10_m=f3m*(h3*ra*P3_cm+l3*(7.5*scmon^2-1.5)*(rm-scmon*ra')');
dx10_s=f3s*(h3*ra*P3_cs+l3*(7.5*scsun^2-1.5)*(rs-scsun*ra')');
dx10=dx10_m+dx10_s; %m

% Contributions to the transverse displacement due to the l(1) term
% IERS10 eq.7.8, diurnal band
    % dt: transverse displacement
dt12_m=-l_11*sin(phi)*f2m*P21_sBm*(sin(phi)*cos(lam-Lm)*n - ...
            cos(2*phi)*sin(lam-Lm)*e);
dt12_s=-l_11*sin(phi)*f2s*P21_sBs*(sin(phi)*cos(lam-Ls)*n - ...
            cos(2*phi)*sin(lam-Ls)*e);
dt12 = dt12_m+dt12_s; %m

dren12=0*r+dt12;

% IERS10 eq.7.9, semidiurnal band
dt13_m=-0.5*l_12*sin(phi)*cos(phi)*f2m*P22_sBm*(cos(2*(lam-Lm))*n + ...
            sin(phi)*sin(2*(lam-Lm))*e);
dt13_s=-0.5*l_12*sin(phi)*cos(phi)*f2s*P22_sBs*(cos(2*(lam-Ls))*n + ...
            sin(phi)*sin(2*(lam-Ls))*e);
dt13 = dt13_m+dt13_s; %m

dren13=0*r+dt13;

% Out of phase contributions from the imaginary parts of h2m and l2m
%IERS10 eq.7.10a, contributions to radial displacement from diurnal tides
dr14_m=-3/4*hI_21*f2m*sin(2*Bm)*sin(2*phi)*sin(lam-Lm);
dr14_s=-3/4*hI_21*f2s*sin(2*Bs)*sin(2*phi)*sin(lam-Ls);
dr14 = dr14_m+dr14_s; %m

% IERS10 eq.7.10b, contributions to transverse displ. from diurnal tides
dt14_m=-3/2*lI_21*f2m*sin(2*Bm)*(cos(2*phi)*sin(lam-Lm)*n + ...
            sin(phi)*cos(lam-Lm)*e);
dt14_s=-3/2*lI_21*f2s*sin(2*Bs)*(cos(2*phi)*sin(lam-Ls)*n + ...
            sin(phi)*cos(lam-Ls)*e);
dt14 = dt14_m+dt14_s; %m

dren14=dr14*r+dt14;

% IERS10 eq.7.11a,contributions to radial displacement from semidiurnal tides
dr15_m=-3/4*hI_22*f2m*(cos(Bm))^2*(cos(phi))^2*sin(2*(lam-Lm));
dr15_s=-3/4*hI_22*f2s*(cos(Bs))^2*(cos(phi))^2*sin(2*(lam-Ls));
dr15 = dr15_m+dr15_s; %m

% IERS10 eq.7.11b, contributions to transverse displ. from semidiurnal tides
dt15_m=3/4*lI_22*f2m*(cos(Bm))^2*(sin(2*phi)*sin(2*(lam-Lm))*n - ...
            2*cos(phi)*cos(2*(lam-Lm))*e);
dt15_s=3/4*lI_22*f2s*(cos(Bs))^2*(sin(2*phi)*sin(2*(lam-Ls))*n - ...
            2*cos(phi)*cos(2*(lam-Ls))*e);
dt15=dt15_m+dt15_s; %m

dren15=dr15*r+dt15;

%%
% Step 2
% Correction due to the frequency variation of the Love and Shida numbers

%%

% Table 7.3a
% DIURNAL tides
% Doodson nr           dRf_ip dRf_op dTf_ip dTf_op

diu =     [1 3 5 6 5 5 -0.08  0.00 -0.01  0.01;     %Q1
           1 4 5 5 4 5 -0.10  0.00  0.00  0.00;
           1 4 5 5 5 5 -0.51  0.00 -0.02  0.03;     %O1
           1 5 5 6 5 5  0.06  0.00  0.00  0.00;     %NO1
           1 6 2 5 5 6 -0.06  0.00  0.00  0.00;     %pi1
           1 6 3 5 5 5 -1.23 -0.07  0.06  0.01;     %P1
           1 6 5 5 4 5 -0.22  0.01  0.01  0.00;     
           1 6 5 5 5 5 12.00 -0.78 -0.67 -0.03;     %K1
           1 6 5 5 6 5  1.73 -0.12 -0.10  0.00;
           1 6 6 5 5 4 -0.50 -0.01  0.03  0.00;     %psi1
           1 6 7 5 5 5 -0.11 -0.11  0.01  0.00];     %phi1
     
Rfid = diu(:,7);
Rfod = diu(:,8);
Tfid = diu(:,9);
Tfod = diu(:,10);


[tau,s,h,p,zns,ps] = doodarg(mjd,leap); %deg
thetafd = diu(:,1)*tau + (diu(:,2)-5)*s + (diu(:,3)-5)*h + ...
            (diu(:,4)-5)*p + (diu(:,5)-5)*zns + (diu(:,6)-5)*ps; %deg
thetafd = thetafd/180*pi; %rad

 dr16 = 0.001* (Rfid.*sin(thetafd+lam) + Rfod.*cos(thetafd+lam))*sin(2*phi);
 dt16 = 0.001*((Tfid.*cos(thetafd+lam) - Tfod.*sin(thetafd+lam))*sin(phi)*e + ...
               (Tfid.*sin(thetafd+lam) + Tfod.*cos(thetafd+lam))*cos(2*phi)*n);    % m

           
dren16=dr16*r+dt16;                      %m
dren16=sum(dren16);
 


%% Eq.17 - long-period band

% Table 7.3b

long = [5 5 5 6 5  0.47  0.16  0.23  0.07
        5 7 5 5 5 -0.20 -0.11 -0.12 -0.05
        6 5 4 5 5 -0.11 -0.09 -0.08 -0.04
        7 5 5 5 5 -0.13 -0.15 -0.11 -0.07
        7 5 5 6 5 -0.05 -0.06 -0.05 -0.03];

Rfil = long(:,6);
Rfol = long(:,7);
Tfil = long(:,8);
Tfol = long(:,9);
   
thetafl = (long(:,1)-5)*s + (long(:,2)-5)*h + (long(:,3)-5)*p + ...
            (long(:,4)-5)*zns + (long(:,5)-5)*ps; %deg
thetafl = thetafl/180*pi; %rad

dr17 = 0.001*(P2_sphi*(Rfil.*cos(thetafl)+Rfol.*sin(thetafl)));        %m
dt17 = 0.001*((Tfil.*cos(thetafl)+Tfol.*sin(thetafl))*sin(2*phi)*n);   %m

dren17=dr17*r+dt17;                 %m
dren17=sum(dren17);


%%

dren12_17=[dren12;dren13;dren14;dren15;dren16;dren17];

% geocentric latitude --> geodetic latitude
%phigd = phigc2phigd(phi,ra);
phigd = cart2phigd(ant);
phi_n(1:6,1)=phigd;
lam_n(1:6,1)=lam;
dx12_17=ren2xyz(dren12_17,phi_n,lam_n);

cts=dx9+dx10+sum(dx12_17);


