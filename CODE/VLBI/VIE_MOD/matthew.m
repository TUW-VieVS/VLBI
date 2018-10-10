% ************************************************************************
%   Description:
%   Solid Earth Tides corrections; displacement caused by lunar and solar
%   gravitational attraction. Translated from Occam matthew.f.
% 
%   Reference: 
%   IERS Conventions, Chapter 7.1.2.
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
%      pdxyz_h             partials for Love numbers in TRF [cm]
%      pdxyz_l             partials for Love numbers in TRF [cm]
% 
%   External calls: 	
%      global rearthm massrSE massrME
%      doodarg.m, cart2phigd., ren2xyz.m, xyz2ren.m
%
%   Coded for VieVS: 
%   10 Oct 2008 by Hana Spicakova
%
%   Revision: 
%   29 Sep 2009 by Lucia Plank: change phi to geodetic latitude
%   25 May 2010 by Hana Spicakova: partials for Love and Shida numbers are
%       transformed into XYZ system (geocentric, terrestrial)
%   12 Nov 2010 by Hana Spicakova: big update concerning frequency
%       dependent Love/Shida numbers. Since now they are computed according
%       to Chapter 7 IERS03 Eq.5b, Table 7.3 and Eq. 8. Values in Tables 7.5a,b
%       are not used anymore
%   12 Nov 2010 by Hana Spicakova: partial derivatives wrt FCN period are
%       added (IERS03 Ch6 Eq.6)
%   11 Apr 2011 by Matthias Madzak: Changed some code for better
%       performance. eg clear i --> i=[];
%   
% *************************************************************************
%
function [cts,pdxyz_h,pdxyz_l,pdxyz_FCN] = matthew(mjd,leap,t2c,ant,xyz_moon,xyz_sun)
  
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
% IERS03 eq.9
% Displacement due to degree-2 tides, with nominal values for h2 and l2
dx9_m=f2m*(h2*ra*P2_cm+3*l2*scmon*(rm-scmon*ra')'); %m
dx9_s=f2s*(h2*ra*P2_cs+3*l2*scsun*(rs-scsun*ra')'); %m
dx9=dx9_m+dx9_s; %m

% IERS03 eq.10
% Displacement due to degree-3 tides
dx10_m=f3m*(h3*ra*P3_cm+l3*(7.5*scmon^2-1.5)*(rm-scmon*ra')');
dx10_s=f3s*(h3*ra*P3_cs+l3*(7.5*scsun^2-1.5)*(rs-scsun*ra')');
dx10=dx10_m+dx10_s; %m

% Contributions to the transverse displacement due to the l(1) term
% IERS03 eq.12, diurnal band
    % dt: transverse displacement
dt12_m=-l_11*sin(phi)*f2m*P21_sBm*(sin(phi)*cos(lam-Lm)*n - ...
            cos(2*phi)*sin(lam-Lm)*e);
dt12_s=-l_11*sin(phi)*f2s*P21_sBs*(sin(phi)*cos(lam-Ls)*n - ...
            cos(2*phi)*sin(lam-Ls)*e);
dt12 = dt12_m+dt12_s; %m

dren12=0*r+dt12;

% IERS03 eq.13, semidiurnal band
dt13_m=-0.5*l_12*sin(phi)*cos(phi)*f2m*P22_sBm*(cos(2*(lam-Lm))*n + ...
            sin(phi)*sin(2*(lam-Lm))*e);
dt13_s=-0.5*l_12*sin(phi)*cos(phi)*f2s*P22_sBs*(cos(2*(lam-Ls))*n + ...
            sin(phi)*sin(2*(lam-Ls))*e);
dt13 = dt13_m+dt13_s; %m

dren13=0*r+dt13;

% Out of phase contributions from the imaginary parts of h2m and l2m
%IERS03 eq.14a, contributions to radial displacement from diurnal tides
dr14_m=-3/4*hI_21*f2m*sin(2*Bm)*sin(2*phi)*sin(lam-Lm);
dr14_s=-3/4*hI_21*f2s*sin(2*Bs)*sin(2*phi)*sin(lam-Ls);
dr14 = dr14_m+dr14_s; %m

% IERS03 eq.14b, contributions to transverse displ. from diurnal tides
dt14_m=-3/2*lI_21*f2m*sin(2*Bm)*(cos(2*phi)*sin(lam-Lm)*n + ...
            sin(phi)*cos(lam-Lm)*e);
dt14_s=-3/2*lI_21*f2s*sin(2*Bs)*(cos(2*phi)*sin(lam-Ls)*n + ...
            sin(phi)*cos(lam-Ls)*e);
dt14 = dt14_m+dt14_s; %m

dren14=dr14*r+dt14;

% IERS03 eq.15a,contributions to radial displacement from semidiurnal tides
dr15_m=-3/4*hI_22*f2m*(cos(Bm))^2*(cos(phi))^2*sin(2*(lam-Lm));
dr15_s=-3/4*hI_22*f2s*(cos(Bs))^2*(cos(phi))^2*sin(2*(lam-Ls));
dr15 = dr15_m+dr15_s; %m

% IERS03 eq.15b, contributions to transverse displ. from semidiurnal tides
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

% Love numbers from the Resonance formulae Ch6, Eq.6
% DIURNAL
% Doodson nr           omega [cpsd]     Hf [mm]
omHf_diu =[1 3 5 6 5 5   0.8908035004645   -50.21
           1 4 5 5 4 5   0.9268492975228   -49.44
           1 4 5 5 5 5   0.9269959896266  -262.25
           1 5 5 6 5 5   0.9638056984426    20.62
           1 6 2 5 5 6   0.9918070331239    -7.16
           1 6 3 5 5 5   0.9945373308561  -122.35
           1 6 5 5 4 5   0.9998514955009    -7.31
           1 6 5 5 5 5   0.9999981876047   369.14
           1 6 5 5 6 5   1.0001448797086    49.97
           1 6 6 5 5 4   1.0027284853370     2.94
           1 6 7 5 5 5   1.0054590443534     5.26
           1 2 5 7 5 5   0.8546110119672    -6.65
           1 2 7 5 5 5   0.8594546483971    -8.02
           1 3 5 6 4 5   0.8906568083607    -9.46
           1 3 7 4 5 5   0.8956471375592    -9.54
           1 4 7 5 5 5   0.9324568463752    -0.10
           1 5 3 6 5 5   0.9583448416940     1.93
           1 5 5 4 5 5   0.9631884787887     7.41
           1 5 5 6 6 5   0.9639523905464     4.13
           1 5 7 4 5 5   0.9686493355373     3.94
           1 6 3 5 4 5   0.9943906387523     1.38
           1 6 4 5 5 4   0.9972676292532     1.02
           1 6 4 5 5 6   0.9972678898725     2.94
           1 6 5 5 7 5   1.0002915718124    -1.07
           1 6 6 5 5 6   1.0027287459563    -0.04
           1 6 6 5 6 4   1.0028750637523     0.08
           1 6 7 3 5 5   1.0048418240346     0.18
           1 7 3 6 5 5   1.0313470396721     3.94
           1 7 5 4 5 5   1.0361906767668    20.62
           1 8 5 5 5 5   1.0730003855829    11.29
           1 8 5 5 6 5   1.0731470776867     7.23];

         
         
% IERS03 Chap 7; Table 7.3
% Parameters in the Resonance Formulae 
%  h_0                   h_2                 l_0                   l_1                   l_2                   l'
prf = ...
[ 0.60671-0.2420e-2i    -0.615e-3-0.122e-4i  0.84963e-1-0.7395e-3i  0.121e-2+0.136e-6i   0.19334e-3-0.3819e-5i -0.221e-3-0.474e-7
 -0.15777e-2-0.7630e-4i  0.160e-5+0.116e-6i -0.22107e-3-0.9646e-5i -0.316e-5-0.166e-6i  -0.50331e-6-0.1639e-7i  0.576e-6+0.303e-7 
  0.18053e-3-0.6292e-5i  0.201e-6+0.279e-8i -0.54710e-5-0.2990e-6i  0.272e-6-0.858e-8i  -0.66460e-8+0.5076e-9i  0.128e-6-0.378e-8
 -0.18616e-5+0.1379e-6i -0.329e-7-0.217e-8i -0.29904e-7-0.7717e-8i -0.545e-8+0.827e-11i  0.10372e-7+0.7511e-9i -0.655e-8-0.291e-9];

om_CW =  -0.0026010-0.0001361i;
om_FCN = 1.0023181+0.000025i;
om_FICN = 0.999026+0.000780i;       


%IERS03 Chap 6 Eq.6
for j=1:size(prf,2)
    x = prf(1,j)+ prf(2,j)./(omHf_diu(:,7)-om_CW) + prf(3,j)./(omHf_diu(:,7)-om_FCN) + prf(4,j)./(omHf_diu(:,7)-om_FICN);
    if j==1
        hfd_0=x;
    elseif j==2
        hfd_2=x;
    elseif j==3
        lfd_0=x;
    elseif j==4
        lfd_1=x;
    elseif j==5
        lfd_2=x;
    elseif j==6
        lfd_c=x;
    end
    % OLD: clear x
    x=[];
end

% OLD: clear i
i=[];

dhf_diu = hfd_0 - (h2+1i*hI_21);
dlf_diu = lfd_0 - (l2+1i*lI_21);

[tau,s,h,p,zns,ps] = doodarg(mjd,leap); %deg
thetafd = omHf_diu(:,1)*tau  + (omHf_diu(:,2)-5)*s   + (omHf_diu(:,3)-5)*h + ...
         (omHf_diu(:,4)-5)*p + (omHf_diu(:,5)-5)*zns + (omHf_diu(:,6)-5)*ps; %deg
thetafd = thetafd/180*pi; %rad

Hfd=omHf_diu(:,8);
 dr16 = 0.001* (real(dhf_diu).*(Hfd*-3/2*sqrt(5/(24*pi))).*sin(thetafd+lam) ...
              + imag(dhf_diu).*(Hfd*-3/2*sqrt(5/(24*pi))).*cos(thetafd+lam))*sin(2*phi);
 
 dt16 = 0.001*((real(dlf_diu).*(Hfd*-3*sqrt(5/(24*pi))).*cos(thetafd+lam) - ...
                imag(dlf_diu).*(Hfd*-3*sqrt(5/(24*pi))).*sin(thetafd+lam))*sin(phi)*e ...
              +(real(dlf_diu).*(Hfd*-3*sqrt(5/(24*pi))).*sin(thetafd+lam) + ...
                imag(dlf_diu).*(Hfd*-3*sqrt(5/(24*pi))).*cos(thetafd+lam))*cos(2*phi)*n);    % m

dren16=dr16*r+dt16;                      %m
dren16=sum(dren16);


%% Eq.17 - long-period band

% LONG PERIOD       
% Doodson nr           omega [cpsd]        Hf [mm]
omHf_long = [5 5 5 6 5  0.000146930783242   27.9
             5 7 5 5 5  0.005461038251366  -30.9
             6 5 4 5 5  0.036192841530055  -35.2
             7 5 5 5 5  0.073002659380692  -66.7
             7 5 5 6 5  0.073148925318761  -27.6];
% Hf is taken from Cartwright D.E. and Tayler R.J. (1971), New computations
% of the Tide-Generating Potential, Geophys.J.R.astr.Soc, 23. pp. 45-74


% IERS03 Chap 7; Eq.8
P = 200/3600/24/366*365; %reference period: 200s [s --> sid.days]
fm = 1/P;       
f = omHf_long(:,6);

% OLD: clear i
i=[];
hfl_0 = 0.5998 - 9.96e-4*(1/(tan(0.15*pi/2))*(1-(fm./f).^0.15) + 1i*(fm./f).^0.15);
lfl_0 = 0.0831 - 3.01e-4*(1/(tan(0.15*pi/2))*(1-(fm./f).^0.15) + 1i*(fm./f).^0.15);
         

long = [5 5 5 6 5  0.47  0.16  0.23  0.07
        5 7 5 5 5 -0.2  -0.11 -0.12 -0.05
        6 5 4 5 5 -0.11 -0.09 -0.08 -0.04
        7 5 5 5 5 -0.13 -0.15 -0.11 -0.07
        7 5 5 6 5 -0.05 -0.06 -0.05 -0.03];
    
Hfl=omHf_long(:,7);

Rfil = sqrt(5/4/pi)*Hfl.*(real(hfl_0)-0.6078);
Rfol = sqrt(5/4/pi)*Hfl.*(-imag(hfl_0)); %out-of-phase contributions from the zonal tides have no closed expression in the time domain
Tfil = 1.5*sqrt(5/4/pi)*Hfl.*(real(lfl_0)-0.0847);
Tfol = 1.5*sqrt(5/4/pi)*Hfl.*(-imag(lfl_0));
  
   
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

%% PARTIALS FOR LOVE AND SHIDA NUMBERS

% Eq. 9
% h^(0)
pdx9_m_h0 = f2m*ra*P2_cm;
pdx9_s_h0 = f2s*ra*P2_cs;
pdx_h0 = pdx9_m_h0+pdx9_s_h0;

% h^(2) - latitude dependence
pdx9_m_h2 = f2m*ra*P2_cm*P2_sphi;
pdx9_s_h2 = f2s*ra*P2_cs*P2_sphi;
pdx_h2 = pdx9_m_h2+pdx9_s_h2;

% l^(0)
pdx9_m_l0 = f2m*3*scmon*(rm-scmon*ra')'; %m
pdx9_s_l0 = f2s*3*scsun*(rs-scsun*ra')'; %m
pdx_l0 = pdx9_m_l0+pdx9_s_l0; %m

% l^(2) - latitude dependence
pdx9_m_l2 = f2m*3*scmon*(rm-scmon*ra')'*P2_sphi; %m
pdx9_s_l2 = f2s*3*scsun*(rs-scsun*ra')'*P2_sphi; %m
pdx_l2 = pdx9_m_l2+pdx9_s_l2; %m

% Eq. 10
% h3
pdx10_m_h3=f3m*ra*P3_cm;
pdx10_s_h3=f3s*ra*P3_cs;
pdx_h3=pdx10_m_h3+pdx10_s_h3; %m

% l3
pdx10_m_l3=f3m*(7.5*scmon^2-1.5)*(rm-scmon*ra')';
pdx10_s_l3=f3s*(7.5*scsun^2-1.5)*(rs-scsun*ra')';
pdx_l3=pdx10_m_l3+pdx10_s_l3; %m


% Eq.12, diurnal band
% l^(1) 
pdt12_m=-sin(phi)*f2m*P21_sBm*(sin(phi)*cos(lam-Lm)*n ...
                            -cos(2*phi)*sin(lam-Lm)*e);
pdt12_s=-sin(phi)*f2s*P21_sBs*(sin(phi)*cos(lam-Ls)*n ...
                            -cos(2*phi)*sin(lam-Ls)*e);
pdt12_l_11 = pdt12_m+pdt12_s; %m



% Eq.13, semidiurnal band
% l^(1)
pdt13_m=-0.5*sin(phi)*cos(phi)*f2m*P22_sBm*(cos(2*(lam-Lm))*n + ...
                                   sin(phi)*sin(2*(lam-Lm))*e);
pdt13_s=-0.5*sin(phi)*cos(phi)*f2s*P22_sBs*(cos(2*(lam-Ls))*n + ...
                                   sin(phi)*sin(2*(lam-Ls))*e);
dt13_l_12 = pdt13_m+pdt13_s; %m


% Eq.14a, out of phase contributions
% hI_21
pdr14_m=-3/4*f2m*sin(2*Bm)*sin(2*phi)*sin(lam-Lm)*r;
pdr14_s=-3/4*f2s*sin(2*Bs)*sin(2*phi)*sin(lam-Ls)*r;
dr_hI_21 = pdr14_m+pdr14_s; %m

% Eq.14b, contributions to transverse displ. from diurnal tides
% lI_21
pdt14_m=-3/2*f2m*sin(2*Bm)*(cos(2*phi)*sin(lam-Lm)*n ...
                           +sin(phi)  *cos(lam-Lm)*e);
pdt14_s=-3/2*f2s*sin(2*Bs)*(cos(2*phi)*sin(lam-Ls)*n ...
                           +sin(phi)  *cos(lam-Ls)*e);
pdt14_lI_21 = pdt14_m+pdt14_s; %m



% Eq.15a, contributions to radial displacement from semidiurnal tides
% hI_22
pdr15_m=-3/4*f2m*(cos(Bm))^2*(cos(phi))^2*sin(2*(lam-Lm))*r;
pdr15_s=-3/4*f2s*(cos(Bs))^2*(cos(phi))^2*sin(2*(lam-Ls))*r;
dr15_hI_22 = pdr15_m+pdr15_s; %m

% Eq.15b, contributions to transverse displ. from semidiurnal tides
% lI_22
pdt15_m=3/4*f2m*(cos(Bm))^2*(sin(2*phi)*sin(2*(lam-Lm))*n ...
                            -2*cos(phi)*cos(2*(lam-Lm))*e);
pdt15_s=3/4*f2s*(cos(Bs))^2*(sin(2*phi)*sin(2*(lam-Ls))*n...
                            -2*cos(phi)*cos(2*(lam-Ls))*e);

dt15_lI_22=pdt15_m+pdt15_s;   %m


% Eq.16, frequency dependence in diurnal band
pdr_dhfR_diu = 0.001* (Hfd*-3/2*sqrt(5/(24*pi))).*sin(thetafd+lam).*sin(2*phi)*r;
pdr_dhfI_diu = 0.001* (Hfd*-3/2*sqrt(5/(24*pi))).*cos(thetafd+lam).*sin(2*phi)*r;
 
pdt_dlfR_diu = 0.001*(Hfd*-3*sqrt(5/(24*pi)).*(cos(thetafd+lam)*sin(phi))*e ...
                    + Hfd*-3*sqrt(5/(24*pi)).*(sin(thetafd+lam)*cos(2*phi))*n);
              
pdt_dlfI_diu = 0.001*(Hfd*-3*sqrt(5/(24*pi)).*(-sin(thetafd+lam)*sin(phi))*e ...
                   +  Hfd*-3*sqrt(5/(24*pi)).*(cos(thetafd+lam)*cos(2*phi))*n);
              
             
           
% Eq.17, frequency dependence in long period band
pdr_Rfil =  0.001*sqrt(5/4/pi)*Hfl*P2_sphi.*cos(thetafl)*r;          %m
pdr_Rfol = -0.001*sqrt(5/4/pi)*Hfl*P2_sphi.*sin(thetafl)*r;  

pdt_Tfil =  0.001*1.5*sqrt(5/4/pi)*Hfl.*cos(thetafl)*sin(2*phi)*n;     %m
pdt_Tfol = -0.001*1.5*sqrt(5/4/pi)*Hfl.*sin(thetafl)*sin(2*phi)*n;    %m


%% Partials for FCN period

ap=(sin(thetafd+lam).*Hfd);
bp=(omHf_diu(:,7)-real(om_FCN)).^2;
cp=ap./bp;
dp=(cos(thetafd+lam).*Hfd);
ep=dp./bp;

pren_FCN = 0.001*(-3/2*sqrt(5/24/pi)*sin(2*phi)*real(prf(3,1))* sum(cp)*r + ...
                -3*sqrt(5/24/pi)*sin(phi)*real(prf(3,3))*sum(ep) *e + ...
                -3*sqrt(5/24/pi)*cos(2*phi)*real(prf(3,3))* sum(cp)*n );
            
%% partials in REN system
pdren_h_14_17=[dr_hI_21; dr15_hI_22; pdr_dhfR_diu; pdr_dhfI_diu; pdr_Rfil; pdr_Rfol];
pdren_l_12_17=[pdt12_l_11; dt13_l_12; pdt14_lI_21; dt15_lI_22; pdt_dlfR_diu; pdt_dlfI_diu; pdt_Tfil; pdt_Tfol];


% transformation from REN into XYZ for h
%clear phi_n lam_n
phi_n=[];
lam_n=[];
phi_n(1:size(pdren_h_14_17,1),1)=phigd;
lam_n(1:size(pdren_h_14_17,1),1)=lam;

pdxyz_h_14_17=ren2xyz(pdren_h_14_17,phi_n,lam_n);

% transformation from REN into XYZ for l
%clear phi_n lam_n
phi_n=[];
lam_n=[];
phi_n(1:size(pdren_l_12_17,1),1)=phigd;
lam_n(1:size(pdren_l_12_17,1),1)=lam;

pdxyz_l_12_17=ren2xyz(pdren_l_12_17,phi_n,lam_n);
%clear phi_n lam_n
%phi_n=[];
%lam_n=[];

% transformation from REN to XYZ
pxyz_FCN=ren2xyz(pren_FCN,phigd,lam); 



% all 77 partials for h in XYZ system (geocentric, terrestrial)
pdxyz_h=[pdx_h0;pdx_h2;pdx_h3;pdxyz_h_14_17]*100; % [cm]

% all 79 partials for l in XYZ system (geocentric, terrestrial)
pdxyz_l=[pdx_l0;pdx_l2;pdx_l3;pdxyz_l_12_17]*100; % [cm]
 
% resonance parameters hRS, lRS, om_FCN
pdxyz_FCN=pxyz_FCN*100; %[cm] 

% pdxyz_h(1): h^(0)
% pdxyz_h(2): h^(2)
% pdxyz_h(3): h3
% pdxyz_h(4): hI_21
% pdxyz_h(5): hI_22
% pdxyz_h(6:36): Rfid
% pdxyz_h(37:67): Rfod
% pdxyz_h(68:72): Rfil
% pdxyz_h(73:77): Rfol

% pdxyz_l(1): l^(0)
% pdxyz_l(2): l^(2)
% pdxyz_l(3): l3
% pdxyz_l(4): l^(1), diurnal
% pdxyz_l(5): l^(1), semidiurnal
% pdxyz_l(6): lI_21
% pdxyz_l(7): lI_22
% pdxyz_l(8:38): Tfid
% pdxyz_l(39:69): Tfod
% pdxyz_l(70:74): Tfil
% pdxyz_l(75:79): Tfol


