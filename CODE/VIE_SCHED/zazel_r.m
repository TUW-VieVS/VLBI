% Purpose  
%   Calculate [az,el] and [ha,dc] with a rigorous model.        
% History  
%   2010-04-06   Jing SUN   Created
%


function [az, el, ha, dc] = zazel_r(mjd, xtrs, ytrs, lon, lat, ra, de, num, ephem)

% 2*pi
twopi = 2 * pi;
% nominal earth rotation velocity [rad/sec]
omega = 7.292115e-5;
% light velocity in m/s 
c = 299792458;    

% trs2crs
[t2c] = ctrs2crs_sched(mjd);

% unit source vector barycentrum-source
SID = sin(de);             
COD = cos(de);             
SIR = sin(ra);
COR = cos(ra);
rq(1) = COD.* COR;   % cos(de) cos(ra)
rq(2) = COD.* SIR;   % cos(de) sin(ra)
rq(3) = SID;         %     sin(de) 
rqu = rq/norm(rq);   % unit source vector barycentrum-source

%
vearth = ephem.earth(num).vbar;
v1 = [-omega*ytrs; omega*xtrs; 0]; % [TRS]
v1 = t2c*v1;                       % [CRS]
k1a = rqu + (vearth+v1)'/c - rqu*((rqu*(vearth+v1))')/c;
% k1a = rqu;
 
% source in TRS (c2t = t2c')
rq_trs = t2c' * k1a';
rq = rq_trs/norm(rq_trs);

% source in local system 
my  = [1,0,0;0,-1,0;0,0,1];  % mirror matrix along y-axis
g2l = my * rotm((pi/2-lat),2) * rotm(lon,3);
lq  = g2l * rq;
   
% zenith distance
zd = acos(lq(3));
% elevation
el = pi/2 - zd;
   
% azimuth(0--2*pi)
saz = atan2(lq(2),lq(1));   % south azimuth
if (lq(2) < 0) 
    saz =  twopi + saz;
end
az = saz + pi;
az = mod(az, twopi);        % north azimuth


% ha (-pi ~ pi)
gmst = tgmst(mjd);
ha = gmst + lon - ra;
while (ha > pi)
    ha = ha - twopi;
end
while (ha < -pi)
    ha = ha + twopi;
end
% dc
dc = de;


