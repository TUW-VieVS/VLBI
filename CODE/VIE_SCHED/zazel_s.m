% Purpose  
%   Calculate [az,el] and [ha,dc] with a simple model. Only the earth 
%   rotation angle as a function of time is taken into account.        
% History  
%   2010-04-06   Jing SUN   Created
%   2016-05-10	Matthias Schartner: improve speed
%   2016-06-14  Matthias Schartner: vectorised funktion
%   2017-06-12  Matthias Schartner: bug-fix for HADC mounts

function [az, el, ha, dc] = zazel_s(mjd, lon, lat, ra, de)
n = length(mjd);

% 2*pi
twopi = 2 * pi;

% earth rotation
tu   = mjd - 51544.5;             % days since fundamental epoch
frac = mjd - floor(mjd) + 0.5;    % Julian day fraction
fac  = 0.00273781191135448;
era  = twopi * (frac + 0.7790572732640 + fac * tu);
era  = mod(era, twopi);           % [rad]

% source vector CRF
sid = sin(de);             
cod = cos(de);             
sir = sin(ra);
cor = cos(ra);
rq = zeros(3,1);
rq(1) = cod * cor;   % cos(de) cos(ra)
rq(2) = cod * sir;   % cos(de) sin(ra)
rq(3) = sid;         %     sin(de) 



% rotation matrix for rotation around z-axis 
caEra = cos(-era);
siEra = sin(-era);

t2c = zeros(n*3,3);
t2c(1:3:end,1,:) = caEra;
t2c(1:3:end,2,:) = -siEra;
t2c(2:3:end,1,:) = siEra;
t2c(2:3:end,2,:) = caEra;
t2c(3:3:end,3,:) = 1;

% source in TRS (c2t = t2c')
rq_trs = t2c*rq;
rq = reshape(rq_trs,3,n);
% rqlenth = sqrt(rq_trs(1,:).^2+rq_trs(2,:).^2+rq_trs(3,:).^2);
% rq     = rq_trs/norm(rq_trs);

% source in local system 
my  = [1,0,0;0,-1,0;0,0,1];  % mirror matrix along y-axis

coLat = cos(pi/2-lat);
siLat = sin(pi/2-lat);
coLon = cos(lon);
siLon = sin(lon);


% rotation matrix for rotation around y-axis and z-axis
% my = eye(3);
% my(2,2)=-1;
% 
% ry = zeros(n*3,3);
% ry(1:3:end,1,:) = coLat;
% ry(1:3:end,3,:) = -siLat;
% ry(2:3:end,2,:) = 1;
% ry(3:3:end,1,:) = siLat;
% ry(3:3:end,3,:) = coLat;
% 
% rz = zeros(n*3,3);
% rz(1:3:end,1,:) = coLon;
% rz(1:3:end,2,:) = siLon;
% rz(2:3:end,1,:) = -siLon;
% rz(2:3:end,2,:) = coLon;
% rz(3:3:end,3,:) = 1;

lq = zeros(3,n);
for i = 1:n
    g2l = [coLat(i) 0 -siLat(i);0 -1 0;siLat(i) 0 coLat(i)]*[coLon(i) siLon(i) 0; -siLon(i) coLon(i) 0;0 0 1];
%     g2l = my*ry(3*i-2:3*i,:)*rz(3*i-2:3*i,:);
    lq(:,i) = g2l*rq(:,i);
end
% g2l = my*[coLat,0,-siLat ;0,1,0;siLat,0,coLat]*[coLon,siLon,0;-siLon,coLon,0;0,0,1];

% zenith distance
zd = acos(lq(3,:));
% elevation
el = pi/2-zd;    

% azimuth(0~2*pi)
saz = atan2(lq(2,:),lq(1,:));   % south azimuth
bool = lq(2,:)<0;
saz(bool)=saz(bool)+twopi;

az = saz + pi;
az = mod(az, twopi);        % north azimuth

       
% ha (-pi ~ pi)
gmst = tgmst(mjd);
ha = gmst + lon - ra;

bool = ha>pi;
while any(bool)
    ha(bool) = ha(bool) - twopi;
    bool = ha>pi;
end
bool = ha<-pi;
while any(bool)
    ha(bool) = ha(bool) + twopi;
    bool = ha<-pi;
end
% dc
dc = de;


