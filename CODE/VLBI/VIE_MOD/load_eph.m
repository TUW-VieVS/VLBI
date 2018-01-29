% ************************************************************************
%   Description:
%   Calculates JPL Ephemerides for VieVs.
%   Comments: ephemerides give barycentric positions in kilometers, output
%   unit is m and m/s; barycenter is solar system barycenter.
% 
%   Reference: 
%   Following the JPL guidelines & the idea of OCCAM 6.2 by O. Titov 
%
%   Input:										
%      tt    (n,1)        time vector Terrestrial time in days MJD (n,1)
%      jplnum (string)    e.g. 'jpl_405' / 'jpl_421'
%                
%   Output:
%      ephem              structure array with barycentric coordinates and
%                         velocities (if selected, also geocentric). Units
%                         are[m, m/s]. Stored at ../EPHEM/ephem.
%                         GM for the planets [m^3/s^2]
% 
%   External calls: 	
%      tt2tdb.m, eph_planet.m
%
%   Loaded data:
%      ../EPHEM/jpl_num.mat
%
%   Coded for VieVS: 
%   18 Aug 2009 by Lucia Plank
%
%   Revision: 
%   05 Mar 2010 by Lucia Plank: enable different JPL version, add GM for 
%                               the planets.
%   08 Aug 2016 by A. Girdiuk: function eph_planet.m is added and new ephemeris jpl_430
%   31 Aug 2016 by A. Girdiuk: Minor change: masses of celestial bodies are stored in the ephem structure also for jpl_430 data.
%
% *************************************************************************
function [ephem]=load_eph(tt,jplnum)

% coonversion TT-time to TDB
  tdb=tt2tdb(tt);

% divide date to last midnight and fraction of day & convert to jd
  et(:,1) = floor(tdb);
  et(:,2) = tdb - et(:,1); 
  et(:,1) = et(:,1) + 2400000.5;
    

%list(i)= 0 --> no interpolation for body i
%       = 1 --> position only
%       = 2 --> position and velocity

 list = [ 2 ... %  = 1: Mercury
          2 ... %  = 2: Venus
          2 ... %  = 3: Earth-Moon Barycenter (needed fo 12)
          2 ... %  = 4: Mars
          2 ... %  = 5: Jupiter
          2 ... %  = 6: Saturn
          2 ... %  = 7: Uranus
          2 ... %  = 8: Neptune
          2 ... %  = 9: Pluto
          2 ... %  =10: Geocentric Moon (needed for 12)
          2 ... %  =11: Sun
          2 ];%   =12: Earth barycentric (needed for geocentric)
      
% unit conversion
posu = 1e3;         % km --> m
velu = 1e3/86400;   % km/day --> m/s

% calculate geocentric positions
geo = 1;

% load ephemerides file
load(strcat('../EPHEM/',jplnum,'.mat'))
 
% error if epoch out of range
err = 0;
if sum(et(1,:)) < jpl.start, err = 1; end
if sum(et(end,:)) > jpl.end, err = 1; end
if err ==1
	disp('ERROR: requested time not within ephemeris limits!!!')
end

% header information
start = jpl.start;
recl  = jpl.reclength;

% get number of record
irec = ceil((et(:,1)- start)./recl);

% fraction of record
frac = ((et(:,1)-((irec-1)*recl + start))+et(:,2))./recl;
dt1  = floor(frac);

% number of sets of coefficients (same for all recs)
nset = jpl.nset;
% number of coefficients per set (same for all recs)
ncoeff = jpl.ncoeff;
% get correct sub-interval number for each body  
%temp = nset*frac;
%iset = floor(temp-dt1)+1;
ephem.time = tt;  % time TT
% normalized chebyshev time (-1 < tc < 1)
%tc   = 2 * (mod(temp,1) + dt1)-1;
% initialize 
pc   = zeros(nset(1),1); pc(1)= 1;
%np   = 2; nv  = 3;
vc   = zeros(nset(1),1); vc(2)= 1;
vfac = (nset + nset)/recl;


for t = 1:length(tt)
  % get correct sub-interval number for each body  
  temp = nset*frac(t);
  iset = floor(temp-dt1(t))+1;
  % normalized chebyshev time (-1 < tc < 1)
  tc   = 2 * (mod(temp,1) + dt1(t))-1;
  
  % ++ calculate for all bodies ++
  % get correct sub-interval number for each body  

	ib =1; % Mercury
	[ephem.merc(t).xbar,ephem.merc(t).vbar]=eph_planet(list(ib),ib,jpl.rec(irec(t)).merc(iset(ib)).set,tc,pc,ncoeff,vc,vfac(ib));

	ib =2; % Venus
	[ephem.venu(t).xbar,ephem.venu(t).vbar]=eph_planet(list(ib),ib,jpl.rec(irec(t)).venu(iset(ib)).set,tc,pc,ncoeff,vc,vfac(ib));

	ib =3; % Earth Moon barycenter
	[ephem.emba(t).xbar,ephem.emba(t).vbar]=eph_planet(list(ib),ib,jpl.rec(irec(t)).emba(iset(ib)).set,tc,pc,ncoeff,vc,vfac(ib));
 
	ib =4; % Mars
	[ephem.mars(t).xbar,ephem.mars(t).vbar]=eph_planet(list(ib),ib,jpl.rec(irec(t)).mars(iset(ib)).set,tc,pc,ncoeff,vc,vfac(ib));
 
	ib =5; % Jupiter
	[ephem.jupi(t).xbar,ephem.jupi(t).vbar]=eph_planet(list(ib),ib,jpl.rec(irec(t)).jupi(iset(ib)).set,tc,pc,ncoeff,vc,vfac(ib));

	ib =6; % Saturn
	[ephem.satu(t).xbar,ephem.satu(t).vbar]=eph_planet(list(ib),ib,jpl.rec(irec(t)).satu(iset(ib)).set,tc,pc,ncoeff,vc,vfac(ib));

	ib =7; % Uranus
	[ephem.uran(t).xbar,ephem.uran(t).vbar]=eph_planet(list(ib),ib,jpl.rec(irec(t)).uran(iset(ib)).set,tc,pc,ncoeff,vc,vfac(ib));

	ib =8; % Neptune
	[ephem.nept(t).xbar,ephem.nept(t).vbar]=eph_planet(list(ib),ib,jpl.rec(irec(t)).nept(iset(ib)).set,tc,pc,ncoeff,vc,vfac(ib));

	ib =9; % Pluto
	[ephem.plut(t).xbar,ephem.plut(t).vbar]=eph_planet(list(ib),ib,jpl.rec(irec(t)).plut(iset(ib)).set,tc,pc,ncoeff,vc,vfac(ib));

	ib =10; % Moon
	[ephem.moon(t).xgeo,ephem.moon(t).vgeo]=eph_planet(list(ib),ib,jpl.rec(irec(t)).moon(iset(ib)).set,tc,pc,ncoeff,vc,vfac(ib));

	ib =11; % Sun
	[ephem.sun(t).xbar,ephem.sun(t).vbar]  =eph_planet(list(ib),ib,jpl.rec(irec(t)).sun(iset(ib)).set, tc,pc,ncoeff,vc,vfac(ib));

	ib =12; % Earth barycentric
    if list(ib)~=0
    	emrat = jpl.emrat;
		ephem.earth(t).xbar = ephem.emba(t).xbar - ephem.moon(t).xgeo/(1+emrat);
		ephem.earth(t).vbar = ephem.emba(t).vbar - ephem.moon(t).vgeo/(1+emrat);
		if list(10)~=0 % barycentric Moon
        	ephem.moon(t).xbar = ephem.moon(t).xgeo + ephem.earth(t).xbar;
			if list(10) == 2 % barycentric Moon velocity
				ephem.moon(t).vbar = ephem.moon(t).vgeo + ephem.earth(t).vbar;
			end
		end
	end
	% geocentric coordinates
	if geo == 1;
		xbar_earth = ephem.earth(t).xbar;
        vbar_earth = ephem.earth(t).vbar;
	    ephem.merc(t).xgeo = ephem.merc(t).xbar - xbar_earth;
		ephem.merc(t).vgeo = ephem.merc(t).vbar - vbar_earth;
		ephem.venu(t).xgeo = ephem.venu(t).xbar - xbar_earth;
		ephem.venu(t).vgeo = ephem.venu(t).vbar - vbar_earth;
	    ephem.mars(t).xgeo = ephem.mars(t).xbar - xbar_earth;
		ephem.mars(t).vgeo = ephem.mars(t).vbar - vbar_earth;
		ephem.jupi(t).xgeo = ephem.jupi(t).xbar - xbar_earth;
		ephem.jupi(t).vgeo = ephem.jupi(t).vbar - vbar_earth;
		ephem.satu(t).xgeo = ephem.satu(t).xbar - xbar_earth;
		ephem.satu(t).vgeo = ephem.satu(t).vbar - vbar_earth;
		ephem.uran(t).xgeo = ephem.uran(t).xbar - xbar_earth;
		ephem.uran(t).vgeo = ephem.uran(t).vbar - vbar_earth;
		ephem.nept(t).xgeo = ephem.nept(t).xbar - xbar_earth;
		ephem.nept(t).vgeo = ephem.nept(t).vbar - vbar_earth;
		ephem.plut(t).xgeo = ephem.plut(t).xbar - xbar_earth;
		ephem.plut(t).vgeo = ephem.plut(t).vbar - vbar_earth;
	    ephem.sun(t).xgeo  =  ephem.sun(t).xbar - xbar_earth;
		ephem.sun(t).vgeo  =  ephem.sun(t).vbar - vbar_earth;
	end
end % time

% conversion factor au^3/day^2 --> m^3/s^2
au = jpl.au*1e3; % km/au ... form jpl header
fau2m = au^3/86400^2;
% gravitational constant GM for the celstial bodies
ephem.gmmerc = jpl.gm(1)*fau2m;
ephem.gmvenu = jpl.gm(2)*fau2m;
ephem.gmmars = jpl.gm(4)*fau2m;
ephem.gmjupi = jpl.gm(5)*fau2m;
ephem.gmsatu = jpl.gm(6)*fau2m;
ephem.gmuran = jpl.gm(7)*fau2m;
ephem.gmnept = jpl.gm(8)*fau2m;
ephem.gmplut = jpl.gm(9)*fau2m;
ephem.gms    = jpl.gm(10)*fau2m;


  
