
%makes the interpolation (spatial with grid_interpol made by Mahdi) and
%temporal of the TEC values, for CODE and igs maps


%function r=interpolation(io,azel)
function r=interpolation(io,latppd,lonppd,h,m,s)

t = h+m/60+s/3600; % time of this scan [h]

if io.lat1<latppd
    TEC=0;  
elseif latppd<io.lat2 
    TEC=0;
elseif io.lon2<lonppd
    TEC=0;
elseif lonppd<io.lon1 
    TEC=0;

else   
    resol = io.interval/3600; % interval of GIM in hours

    for i = 1 : io.number_of_maps-1
        if resol*(i-1) <= t && t <resol*i
            TEC1=grid_interpol(lonppd,latppd,io,i);
            TEC2=grid_interpol(lonppd,latppd,io,i+1);
            
            TEC=(TEC1+(((t- (resol*(i-1)) )*(TEC2-TEC1))/2))*10^(io.exponent);     
        end
    end
 end
   
 r=TEC;