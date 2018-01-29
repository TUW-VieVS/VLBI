
%makes the interpolation (spatial with grid_interpol made by Mahdi) and
%temporal of the TEC values, for CODE and igs maps


%function r=interpolation(io,azel)
function r=interpolation(io,latppd,lonppd,h,m,s)


if io.lat1<latppd
        TEC=0;  else 
        if latppd<io.lat2 
         TEC=0;  else
         if io.lon2<lonppd
        TEC=0; else 
        if lonppd<io.lon1 
         TEC=0;  else
   
     if 0<=h&&h<2
        TEC1=grid_interpol(lonppd,latppd,io,1);
        TEC2=grid_interpol(lonppd,latppd,io,2);
        t=h+m/60+s/3600;
        TEC=(TEC1+(((t-0)*(TEC2-TEC1))/2))*10^(io.exponent);  else
     if 2<=h&&h<4
        TEC1=grid_interpol(lonppd,latppd,io,2);
        TEC2=grid_interpol(lonppd,latppd,io,3);
        t=h+m/60+s/3600;
        TEC=(TEC1+(((t-2)*(TEC2-TEC1))/2))*10^(io.exponent);  else
    if 4<=h&&h<6
        TEC1=grid_interpol(lonppd,latppd,io,3);
        TEC2=grid_interpol(lonppd,latppd,io,4);
        t=h+m/60+s/3600;
        TEC=(TEC1+(((t-4)*(TEC2-TEC1))/2))*10^(io.exponent);    else   
    if 6<=h&&h<8
        TEC1=grid_interpol(lonppd,latppd,io,4);
        TEC2=grid_interpol(lonppd,latppd,io,5);
        t=h+m/60+s/3600;
        TEC=(TEC1+(((t-6)*(TEC2-TEC1))/2))*10^(io.exponent);  else      
    if 8<=h&&h<10
        TEC1=grid_interpol(lonppd,latppd,io,5);
        TEC2=grid_interpol(lonppd,latppd,io,6);
        t=h+m/60+s/3600;
        TEC=(TEC1+(((t-8)*(TEC2-TEC1))/2))*10^(io.exponent);    else     
    if 10<=h&&h<12
        TEC1=grid_interpol(lonppd,latppd,io,6);
        TEC2=grid_interpol(lonppd,latppd,io,7);
        t=h+m/60+s/3600;
        TEC=(TEC1+(((t-10)*(TEC2-TEC1))/2))*10^(io.exponent);   else     
    if 12<=h&&h<14
        TEC1=grid_interpol(lonppd,latppd,io,7);
        TEC2=grid_interpol(lonppd,latppd,io,8);
        t=h+m/60+s/3600;
        TEC=(TEC1+(((t-12)*(TEC2-TEC1))/2))*10^(io.exponent);   else     
    if 14<=h&&h<16
        TEC1=grid_interpol(lonppd,latppd,io,8);
        TEC2=grid_interpol(lonppd,latppd,io,9);
        t=h+m/60+s/3600;
        TEC=(TEC1+(((t-14)*(TEC2-TEC1))/2))*10^(io.exponent);      else  
    if 16<=h&&h<18
        TEC1=grid_interpol(lonppd,latppd,io,9);
        TEC2=grid_interpol(lonppd,latppd,io,10);
        t=h+m/60+s/3600;
        TEC=(TEC1+(((t-16)*(TEC2-TEC1))/2))*10^(io.exponent);   else     
    if 18<=h&&h<20
        TEC1=grid_interpol(lonppd,latppd,io,10);
        TEC2=grid_interpol(lonppd,latppd,io,11);
        t=h+m/60+s/3600;
        TEC=(TEC1+(((t-18)*(TEC2-TEC1))/2))*10^(io.exponent); else       
    if 20<=h&&h<22
        TEC1=grid_interpol(lonppd,latppd,io,11);
        TEC2=grid_interpol(lonppd,latppd,io,12);
        t=h+m/60+s/3600;
        TEC=(TEC1+(((t-20)*(TEC2-TEC1))/2))*10^(io.exponent);    else    
    if 22<=h
        TEC1=grid_interpol(lonppd,latppd,io,12);
        TEC2=grid_interpol(lonppd,latppd,io,12);
        t=h+m/60+s/3600;
        TEC=(TEC1+(((t-22)*(TEC2-TEC1))/2))*10^(io.exponent); 
  
    
   
    
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
 end
    
 r=TEC;