% ************************************************************************
%   Description:
%   Gives the local elevation angle as a function of time, station 
%   coordinates and declination and right ascension of the source. Nutation
%   polar motion are not included. Only the earth rotation angle as a
%   function of time is taken into account. Tests showed that the result is
%   accurate to the 0.1 degree level. The function is called from read_ngs,
%   in order to eliminate observations below the cutoff angle.
%
%   Input:										
%     mjd               time in MJD [d]
%     ant               antenna coordinates (x,y,z), [m]
%     de                declination of the source [rad]
%     ra                rigth ascension of the source [rad]
% 
%   Output:
%     el                local elevation angle [rad]
% 
%   External calls: 	
%   rotm.m, xyz2ell.m     
%       
%   Coded for VieVS: 
%   02 Dec 2009 by Lucia Plank
%
%   Revision: 
%   17 Dec 2009 by Lucia Plank: replace xyz2ellip.m with xyz2ell.m
%
% ************************************************************************
function el = elev(mjd,ant,de,ra,varargin)

 tu   = mjd - 51544.5;             % days since fundamental epoch
 frac = mjd - floor(mjd) + 0.5;    % Julian day fraction
 fac  = 0.00273781191135448;
 era  = 2*pi * (frac + 0.7790572732640 + fac * tu );
 era  = mod(era,2*pi);              % [rad]
 
 % rotation: 
 %     angle: era 
 %     axes: z
 cera = cos(-era);
 sera = sin(-era);
 t2c  = [cera, sera 0; -sera cera 0; 0 0 1];
%  rotm(-era,3);

 % source vector CRF
 
 sid = sin(de);             
 cod = cos(de);             
 sir = sin(ra);
 cor = cos(ra);

 rq(1) = cod * cor;           % cos(de) cos(ra)
 rq(2) = cod * sir;           % cos(de) sin(ra)
 rq(3) = sid;                 %     sin(de) 
 
 % ellipsoidal antenna coordinates       
 if ~isempty(varargin)
    phi = varargin{1};
    lam = varargin{2};
 else
   [phi,lam] =xyz2ell(ant);   %[rad,rad,m]
 end
 % source in TRS (c2t = t2c')
   rq_trs = t2c'*rq';
   rq     = rq_trs/norm(rq_trs);

   % source in local system 
   my  = [1,0,0;0,-1,0;0,0,1];  % mirror matrix along y-axis
   c2 = cos(pi/2-phi);
   s2 = sin(pi/2-phi);
   clam = cos(lam);
   slam = sin(lam);
   
 % rotation: 
 %     angle: pi/2-phi 
 %     axes: y
 % rotation: 
 %     angle: lam 
 %     axes: z
   g2l = my*[c2 0 -s2; 0 1 0; s2 0 c2]*[clam slam 0; -slam clam 0; 0 0 1];
%    rotm((pi/2-phi),2)*rotm(lam,3);
   
   lq = g2l*rq;

   % zenith distance
   zd = acos(lq(3));
   el = pi/2-zd;    
    
    
    