% -------------------------------
% IONEX bivariate interpolation
% -------------------------------
% M. Mahdi Alizadeh
%  03.04.2010
%
%   04 May 2011 by Matthias Madzak: Preallocating variable 'map' for better
%       performance.
% -------------------------------

function E = grid_interpol(lambda,beta,ionex,map_one_h)

% + 04 May 2011 by Matthias Madzak
% preallocate
map=zeros(size(ionex.map));
% - 04 May 2011 by Matthias Madzak

if ionex.lat1>ionex.lat2
    lat1=ionex.lat2;
    lat2=ionex.lat1;
    dlat=-ionex.dlat;
    for i=1:ionex.number_of_maps
        map(:,:,i)=flipud(ionex.map(:,:,i));        
    end
else
    lat1=ionex.lat1;
    lat2=ionex.lat2;
    dlat=ionex.dlat;
    for i=1:ionex.number_of_maps
        map(:,:,i)=ionex.map(:,:,i);
    end
end

if ionex.lon1>ionex.lon2
    lon1=ionex.lon2;
    lon2=ionex.lon1;
    dlon=-ionex.dlon;
    for i=1:ionex.number_of_maps
        map(:,:,i)=fliplr(ionex.map(:,:,i));
    end
else
    lon1=ionex.lon1;
    lon2=ionex.lon2;
    dlon=ionex.dlon;
end


ion_lon = lon1:dlon:lon2;
pos_lambda = find(lambda > ion_lon);
d1 = length(pos_lambda);
lambda_point_00 = ion_lon(d1);

ion_lat = lat1:dlat:lat2;
pos_beta = find(beta > ion_lat);
d2 = length(pos_beta);
beta_point_00 = ion_lat(d2);

p = (lambda - lambda_point_00)/5;
q = (beta - beta_point_00)/2.5;


E_00 = map(d2,d1,map_one_h);
E_10 = map(d2,d1+1,map_one_h);
E_01 = map(d2+1,d1,map_one_h);
E_11 = map(d2+1,d1+1,map_one_h);

E = (1-p)*(1-q)*E_00 + p*(1-q)*E_10 + q*(1-p)*E_01 + p*q*E_11;

