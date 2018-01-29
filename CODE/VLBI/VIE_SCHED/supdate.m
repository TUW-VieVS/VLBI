% Purpose  
%   Update the scheduling results.
% History  
%   2010-04-06   Jing SUN   Created
%   
% Changes
%   2015-07-15, A. Hellerschmied: Comments added.
%   2016-05-24, M. Schartner: bugfix
%   2016-06-21, M. Schartner: now saves az, el, ha, dc at endmjd
%   2016-07-11, M. Schartner: possibility to dont save the observation for check with PARA.MIN_SRCRP

function [staobs, srcobs, endmjdmax, scanuminfo] = supdate(station, source, subcon, staobs, srcobs, scanuminfo, PARA, varargin)

% do you want to update sky coverage?
saveSource = true;
output = true;
if length(varargin)>=1
    if any(strcmp(varargin,'dontSaveSource'))
    saveSource = false;
    end
    if any(strcmp(varargin,'noOutput'))
    output = false;
    end
end

%
stanum = length(station);
srcnum = length(source);

% Calc az el ha dc at endmjd
azEndmjd = zeros(stanum,1);
elEndmjd = zeros(stanum,1);
haEndmjd = zeros(stanum,1);
dcEndmjd = zeros(stanum,1);
% get az el ha dc at startmjd
azStartmjd = zeros(stanum,1);
elStartmjd =zeros(stanum,1);

for iscan = 1 : subcon.nscan
    allsta = [subcon.scan(iscan).sta];
    allsta_endmjd = [allsta.endmjd];
    staid = [allsta.staid];
    llh = [station(staid).llh];
    lon = llh(1:3:end);
    lat = llh(2:3:end);
    ra = source(subcon.scan(iscan).srcid).ra;
    de = source(subcon.scan(iscan).srcid).de;
    [az, el, ha, dc] = zazel_s(allsta_endmjd, lon, lat, ra, de);
    azEndmjd(staid)=az;
    elEndmjd(staid)=el;
    haEndmjd(staid)=ha;
    dcEndmjd(staid)=dc;
    
    azStartmjd(staid)=[allsta.az];
    elStartmjd(staid)=[allsta.el];
end

% compare az al ha dc startmjd with endmjd
diffAz = azStartmjd-azEndmjd;
diffAzRound = round(diffAz/(2*pi));
azEndmjd = azEndmjd+diffAzRound*(2*pi);

% are there problems?
% if max(abs(azStartmjd-azEndmjd))>pi/4
%     error('error with az cable wrap!')
% end
% if max(abs(elStartmjd-elEndmjd))>pi/8
%     error('error with el cable wrap!')
% end



% #### update staobs ####
for iscan = 1 : subcon.nscan
    srcid = subcon.scan(iscan).srcid;
    for ista = 1 : subcon.scan(iscan).nsta 
        staid = subcon.scan(iscan).sta(ista).staid;
        staobs(staid).ifirst = 0;
        staobs(staid).endmjd = subcon.scan(iscan).sta(ista).endmjd;
        staobs(staid).az     = azEndmjd(staid);
        staobs(staid).el     = elEndmjd(staid);
        staobs(staid).ha     = haEndmjd(staid);
        staobs(staid).dc     = dcEndmjd(staid);
        
        [cat] = zcoverage(staobs(staid).az, staobs(staid).el); % Calculate the sky coverage number for the current station
        staobs(staid).ncat = staobs(staid).ncat + 1;
        staobs(staid).cattn(staobs(staid).ncat) = staobs(staid).endmjd;
        staobs(staid).catsn(staobs(staid).ncat) = cat;
        % calc new sky-coverage vector. 
        ncat = 0;
        cattn = [];
        catsn = [];
        for ic = 1 : staobs(staid).ncat
            % PARA.SKYDT = The interval for calculation of sky coverage [min] 
            if (((staobs(staid).endmjd-staobs(staid).cattn(ic))*1440) <= PARA.SKYDT) 
                ncat = ncat + 1;
                cattn(ncat) = staobs(staid).cattn(ic);
                catsn(ncat) = staobs(staid).catsn(ic);
            end
        end
        % update new sky coverage vector
        staobs(staid).ncat = ncat;
        staobs(staid).cattn = cattn;
        staobs(staid).catsn = catsn;
        
        
    end
end

% #### calculate endmjdmax ####
endmjdmax = 0.0;
for ista = 1 : stanum
    if (staobs(ista).endmjd > endmjdmax)
        endmjdmax = staobs(ista).endmjd;
    end
end

% #### check the downtime ####
for ista = 1 : stanum
    for idn = 1 : station(ista).downum
        if ((staobs(ista).endmjd >= station(ista).downstart(idn)) & (staobs(ista).endmjd <= station(ista).downend(idn)))
            staobs(ista).endmjd = endmjdmax;
        end
    end
end 
for ista = 1 : stanum
    for idn = 1 : station(ista).downum
        if ((endmjdmax >= station(ista).downstart(idn)) & (endmjdmax <= station(ista).downend(idn)))
            staobs(ista).endmjd = endmjdmax;
        end
    end
end 

% #### update srcobs ####
for iscan = 1 : subcon.nscan
    srcid = subcon.scan(iscan).srcid;
    if saveSource
        srcobs(srcid).obsmjd = subcon.scan(iscan).startmjd;
        srcobs(srcid).nscan = srcobs(srcid).nscan + 1;
    end
end

% #### print info ####
if (PARA.SCREEN == 1) && output
    for iscan = 1 : subcon.nscan
        scanuminfo = scanuminfo + 1;
        fprintf(PARA.fid_body,'scan %3d ', scanuminfo);
        
        if ~strcmp(source(subcon.scan(iscan).srcid).commoname,'        ')
            fprintf(PARA.fid_body,' %8s ', source(subcon.scan(iscan).srcid).commoname);
        else
            fprintf(PARA.fid_body,' %8s ', source(subcon.scan(iscan).srcid).name);
        end
        
        [year, month, day, hour, minute, second] = tymdhms(subcon.scan(iscan).startmjd);
        [datestr] = tdatestr(year, month, day, hour, minute, second);    
        fprintf(PARA.fid_body,' %s   #sta: %d \n', datestr,subcon.scan(iscan).nsta);
        for ista = 1 : subcon.scan(iscan).nsta
            staid = subcon.scan(iscan).sta(ista).staid;
            waitsec = ceil((subcon.scan(iscan).startmjd-subcon.scan(iscan).sta(ista).startmjd)*86400);
            fprintf(PARA.fid_body,' %s %4dslew-->%4didle-->%4dobs\n', station(staid).name, subcon.scan(iscan).sta(ista).slewtime, waitsec, subcon.scan(iscan).sta(ista).duration);
        end
    end
end

%
loopunit = 1;   %[minute]
if(subcon.nscan == 0)
    for ista = 1 : stanum
        staobs(ista).endmjd = staobs(ista).endmjd + loopunit/1440;
    end
end


