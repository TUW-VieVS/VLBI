% loads vEarth frot JPL Ephemerides at time tt
% 04-07-2016 M. Schartner: created

function [ ephem ] = get_vEarth( tt,jpl )
% coonversion TT-time to TDB
tdb=tt2tdb(tt);

% divide date to last midnight and fraction of day & convert to jd
et(:,1) = floor(tdb);
et(:,2) = tdb - et(:,1); 
et(:,1) = et(:,1) + 2400000.5;
    

%list(i)= 0 --> no interpolation for body i
%       = 1 --> position only
%       = 2 --> position and velocity

 list = [ 2 ... %  = 1: Earth-Moon Barycenter (needed for 3)
          2 ... %  = 2: Geocentric Moon (needed for 3)
          2 ];%   = 3: Earth barycentric (needed for geocentric)
      
% unit conversion
posu = 1e3;         % km --> m
velu = 1e3/86400;   % km/day --> m/s

% calculate geocentric positions
geo = 1;
 
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


for t = 1:length(tt)
    ephem.time(t) = tt(t);  % time TT
    % get correct sub-interval number for each body  
    temp = nset*frac(t);
    iset = floor(temp-dt1(t))+1;
    % normalized chebyshev time (-1 < tc < 1)
    tc   = 2 * (mod(temp,1) + dt1(t))-1;
    % initialize 
    pc   = zeros(nset(1),1); pc(1)= 1;
    np   = 2; nv  = 3;
    vc   = zeros(nset(1),1); vc(2)= 1;
    vfac = (nset + nset)/recl;
    twot = 0;
  
    % ++ calculate for all bodies ++
    ib =1; % Earth Moon barycenter
	% calculate polynomial values
    if tc(ib) ~= pc(2)
        np    = 2; nv = 3;
        pc(2) = tc(ib);
        twot  = tc(ib) + tc(ib);
    end
    if np < ncoeff(ib)
         for i = np+1:ncoeff(ib)
            pc(i) = twot*pc(i-1)- pc(i-2);
         end
         np=ncoeff;
    end
    % interpolate and get position for each component  
    set = jpl.rec(irec(t)).emba(iset(ib)).set;
    ephem.emba(t).xbar = (pc(1:ncoeff(ib))'* set)'*posu;
    % Velocity
    vc(3)= twot + twot;    
    if nv < ncoeff(ib)
        for i=nv+1:ncoeff(ib)     
            vc(i) = twot*vc(i-1) + pc(i-1) + pc(i-1) - vc(i-2);
        end
        nv = ncoeff(ib);
    end
    %  iterpolate and get velocity for each component
    vel = (vc(1:ncoeff(ib))'*set)';
    ephem.emba(t).vbar = vfac(ib) * vel*velu;    

    ib =2; % Moon
    % calculate polynomial values
    if tc(ib) ~= pc(2)
        np    = 2; nv = 3;
        pc(2) = tc(ib);
        twot  = tc(ib) + tc(ib);
    end
    if np < ncoeff(ib)
         for i = np+1:ncoeff(ib)
            pc(i) = twot*pc(i-1)- pc(i-2);
         end
         np=ncoeff;
    end
    % interpolate and get position for each component  
    set = jpl.rec(irec(t)).moon(iset(ib)).set;
    ephem.moon(t).xgeo = (pc(1:ncoeff(ib))'* set)'*posu;
    % Velocity
    vc(3)= twot + twot;    
    if nv < ncoeff(ib)
        for i=nv+1:ncoeff(ib)     
            vc(i) = twot*vc(i-1) + pc(i-1) + pc(i-1) - vc(i-2);
        end
    nv = ncoeff(ib);
    end
    %  iterpolate and get velocity for each component
    vel = (vc(1:ncoeff(ib))'*set)';
    ephem.moon(t).vgeo = vfac(ib) * vel*velu;    
    
    ib =3; % Earth barycentric
    if list(ib)~=0
        emrat = jpl.emrat;
        ephem.earth(t).xbar = ephem.emba(t).xbar - ephem.moon(t).xgeo/(1+emrat);
        ephem.earth(t).vbar = ephem.emba(t).vbar - ephem.moon(t).vgeo/(1+emrat);
        ephem.moon(t).xbar = ephem.moon(t).xgeo + ephem.earth(t).xbar;
        ephem.moon(t).vbar = ephem.moon(t).vgeo + ephem.earth(t).vbar;
    end
end 

end

