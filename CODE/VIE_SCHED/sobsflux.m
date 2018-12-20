% Purpose  
%   Calculate observed flux density of source.
% History  
%   2010-04-06   Jing SUN   created
%   2016-05-10	Matthias Schartner: improve speed slightly
% 	2016-07-06						vectorised function									



function [obsflux] = sobsflux(mjd, blx, bly, blz, ra, de, fluxpara, PARA)

% constant
flcon1 = (pi * pi) / (4.0 * 0.6931471);
flcon2 = pi / (3600.0 * 180.0 * 1000.0);

%
gmst = tgmst(mjd);
cosdec = cos(de);
ha = gmst - ra;

u = blx * sin(ha) + bly * cos(ha);
v = blz*cos(de) + sin(de) .* (-blx * cos(ha) + bly * sin(ha));

n = length(blx);
% compute observed flux density
obsflux = zeros(n,PARA.MAX_BANDNUM);
for iband = 1 : PARA.MAX_BANDNUM  
    if (abs(fluxpara(iband,PARA.MAX_FLUXPARA)-1.0)<1.0d-6)   % type B
        % compute projected baseline
        pbase = sqrt(u.*u + v.*v) / 1000.0;
        boolsave = false(n,1);
        for i = 1 : (PARA.MAX_FLUXPARA/2-1)   % up to eight fluxes and nine baseline lengths may be specified 
            ib = i * 2;
            
            if fluxpara(iband,ib) > 0.0
            bool1 = pbase >= fluxpara(iband,ib-1);
            bool2 = pbase <= fluxpara(iband,ib+1);
            bool = bool1 & bool2 ;
            boolsave = boolsave | bool;
            obsflux(bool,iband) = fluxpara(iband,ib);
            
            if all(boolsave)
                break
            end
            end
        end
    elseif (abs(fluxpara(iband,PARA.MAX_FLUXPARA)-2.0)<1.0d-6)   % type M
        % compute u,v
        u = u / PARA.WAVEL(iband);  % from m to wavelengths
        v = v / PARA.WAVEL(iband);  % from m to wavelengths
        % cpnum (up to three components)
        if (fluxpara(iband,13) > 1.0d-3)   
            cpnum = 3;
        elseif (fluxpara(iband,7) > 1.0d-3)  
            cpnum = 2;
        elseif (fluxpara(iband,1) > 1.0d-3)  
            cpnum = 1;
        end   
        for i = 1 : cpnum    
            im = (i - 1) * 6 + 1; 
            fluxmax = fluxpara(iband,im);               % [Jy]
            majax = fluxpara(iband,(im+1)) * flcon2;    % radians
            ratio = fluxpara(iband,(im+2));             % no units
            pa = fluxpara(iband,(im+3)) * pi / 180.0;   % radians
            ucospa = u * cos(pa);
            usinpa = u * sin(pa);
            vcospa = v * cos(pa);
            vsinpa = v * sin(pa);
            arg1 = (vcospa + usinpa) .* (vcospa + usinpa);
            arg2 = (ratio * (ucospa - vsinpa)) .* (ratio * (ucospa - vsinpa));
            arg = -flcon1 * (arg1 + arg2) * majax * majax;
            fl = fluxmax * exp(arg);
            obsflux(:,iband) = obsflux(:,iband) + fl;
        end       
    end
end   


