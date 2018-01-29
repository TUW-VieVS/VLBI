% Purpose  
%   Calculate the possible subconfiguration (a group of scans).
% History  
%   2010-04-06   Jing SUN   Created
% 
% Changes
%   2015-07-15, A. Hellerschmied: Comments added.
%   2015-10-28, A. Hellerschmied: Comments added.
%   2016-05-10,	Matthias Schartner: improve speed slightly
%                                   added STAR mode
%   2016-06-07, Matthias Schartner: vectorised sobsflux
%   2016-07-11, Matthias Schartner: small changes to improve readability 
%   2016-09-08, Matthias Schartner: bugfix 


function [subcon] = ssubcon_sta(source, srcobs, station, staobs, PARA, varargin)

% optional 6th parameter is used for STAR mode  
if nargin == 6 
    starsource = varargin{1};
else
    starsource = 0;
end

subcon(1).nscan = 0;
srcnum = length(source);
stanum = length(station);
% endmjdmax
endmjdmax = 0.0;
for ista = 1 : stanum
    if (staobs(ista).endmjd > endmjdmax)
        endmjdmax = staobs(ista).endmjd;
    end
end

llh = [station.llh];
lon = llh(1:3:end);
lat = llh(2:3:end);
% check the visibility at stations for all sources
for isrc = 1 : srcnum
    ra = source(isrc).ra;
    de = source(isrc).de;
    srcsta(isrc).downstanum = 0;            % Number of stations where the current source is UP (visible)!
    srcsta(isrc).downstasn(1:stanum) = 0;   % List of stations where the current source is UP (visible)!
    endmjd = [staobs.endmjd];
    
    [az, el, ha, dc] = zazel_s(endmjd, lon, lat, ra, de);
    for ista = 1 : stanum
        [lup] = zlup(endmjd(ista), az(ista), el(ista), ha(ista), dc, ista, station, PARA.MIN_CUTEL);
        if (lup == 1)
            srcsta(isrc).downstanum = srcsta(isrc).downstanum + 1;
            srcsta(isrc).downstasn(srcsta(isrc).downstanum) = ista;
        end
    end  
end

% initialize not possible struct
np = struct;

% ##### check the possible subconfiguration with one source #####
nsub = 1;
for isrc = 1 : srcnum
    np(isrc).notPossible=[];
    np(isrc).bool = [];

    % check the visibility
    if (srcsta(isrc).downstanum < max(2,PARA.MIN_STANUM))
        continue;
    end
    % check if the source is repeated in a short interval 
    if ((endmjdmax-srcobs(isrc).obsmjd) < PARA.MIN_SRCRP/1440)
        continue;
    end 
    %check source projected flux
    ra = source(isrc).ra;
    de = source(isrc).de;
    fluxpara(1:PARA.MAX_BANDNUM,1:PARA.MAX_FLUXPARA) = source(isrc).fluxpara(1:PARA.MAX_BANDNUM,1:PARA.MAX_FLUXPARA);
    
    conf = nchoosek(srcsta(isrc).downstasn(1:srcsta(isrc).downstanum),2);

    staid1 = conf(:,1);
    staid2 = conf(:,2); 
    
    xyz1 = [station(staid1).xyz]';
    xyz2 = [station(staid2).xyz]';
    
    blx = xyz1(:,1)-xyz2(:,1);
    bly = xyz1(:,2)-xyz2(:,2);
    blz = xyz1(:,3)-xyz2(:,3);

    [obsflux] = sobsflux(endmjdmax, blx, bly, blz, ra, de, fluxpara, PARA); 
    bool = obsflux(:,1:PARA.MAX_BANDNUM) >= PARA.MIN_FLUX;
    % which baseline1 can be observed?
    bool = sum(bool,2)==size(bool,2);
    notPossible = [staid1(~bool,:) staid2(~bool,:)];
    np(isrc).notPossible=notPossible;
    np(isrc).bool = bool;
    if source(isrc).star == 0 && (min(obsflux(:)) < PARA.MIN_FLUX)  % Observed flux large enough?
        continue; % Check next source
    end
    
    % check if source and STAR mode matches. If no STAR mode is used, value
    % of starsource and source(isrc).star is always 0.
    if starsource==1 && source(isrc).star == 0   
        continue;
    end
    if starsource==0 && source(isrc).star == 1 
        continue;
    end 
    
    % the possible subconfiguration
    subcon(nsub).nscan = 1;
    subcon(nsub).scan(1).srcid = isrc;
    subcon(nsub).scan(1).startmjd = 0.0;
    subcon(nsub).scan(1).nsta = srcsta(isrc).downstanum;
    for ista = 1 : subcon(nsub).scan(1).nsta
        subcon(nsub).scan(1).sta(ista).staid    = srcsta(isrc).downstasn(ista);
        subcon(nsub).scan(1).sta(ista).slewtime = 0.0;
        subcon(nsub).scan(1).sta(ista).startmjd = 0.0;
        subcon(nsub).scan(1).sta(ista).az       = 0.0;
        subcon(nsub).scan(1).sta(ista).el       = 0.0;
        subcon(nsub).scan(1).sta(ista).ha       = 0.0;
        subcon(nsub).scan(1).sta(ista).dc       = 0.0;
        subcon(nsub).scan(1).sta(ista).duration = 0.0;
        subcon(nsub).scan(1).sta(ista).endmjd   = 0.0;
    end 
    nsub = nsub + 1;
end


if PARA.SUBNETTING == 2
    if starsource == 0
        % ##### check the possible subconfiguration with two sources #####
        for isrc1 = 1 : srcnum-1
            if starsource==0 && source(isrc1).star == 1   % skip source if observation is for STAR mode
                continue;
            end
            downstanum1 = srcsta(isrc1).downstanum;
            downstasn1 = srcsta(isrc1).downstasn;

            % get baselines which are not possible
            notPossible1 = np(isrc1).notPossible;
            bool1 = np(isrc1).bool;
            if all(~bool1) % no baseline possible
                continue
            end

            % check the visibility
            if (downstanum1 < max(2,PARA.MIN_STANUM)) 
                continue;
            end

            if endmjdmax-srcobs(isrc1).obsmjd < PARA.MIN_SRCRP/1440
                continue;
            end 

            for isrc2 = isrc1+1 : srcnum
                if starsource==0 && source(isrc2).star == 1   % skip source if observation is for STAR mode
                    continue;
                end
                downstanum2 = srcsta(isrc2).downstanum;
                downstasn2 = srcsta(isrc2).downstasn;

                % check the visibility
                if (downstanum2 < max(2,PARA.MIN_STANUM))
                    continue;
                end

                % check if the source is repeated in a short interval 
                if endmjdmax-srcobs(isrc2).obsmjd < PARA.MIN_SRCRP/1440
                    continue;
                end 

                % check the angle distance
                ra1 = source(isrc1).ra;
                de1 = source(isrc1).de;
                ra2 = source(isrc2).ra;
                de2 = source(isrc2).de;
                arg = cos(de1) * cos(de2) * cos(ra1 - ra2) + sin(de1) * sin(de2);
                srcd = acos(arg);
                if (srcd < PARA.MIN_SRC2ANG)
                    continue;
                end 

                notPossible2 = np(isrc2).notPossible;
                bool2 = np(isrc2).bool;
                % if no baseline is working
                if all(~bool2) 
                    continue
                end

                % calculate the possible grouping
                [~, pmcb2] = spmcb2(downstanum1, downstasn1, downstanum2, downstasn2);
                pmcb2([pmcb2.downstanum1]==1 | [pmcb2.downstanum2]==1)=[];
                pmcb2([pmcb2.downstanum1]<PARA.MIN_STASCAN | [pmcb2.downstanum2]<PARA.MIN_STASCAN)=[];
                np2 = length(pmcb2);

                % Loop over all groupings:
                for ip = 1 : np2

                    % ##### Source 1 #####
                    sta1 = pmcb2(ip).downstasn1;
                    if length(sta1)<2
                        continue
                    end
                    flag = false;
                    for i = 1:size(notPossible1,1)
                        if all(ismember(notPossible1(i,:),sta1))
                            flag = true;
                            break;
                        end   
                    end

                    if flag
                        continue
                    end

                    % ##### Source 2 #####
                    sta2 = pmcb2(ip).downstasn2;
                    if length(sta2)<2
                        continue
                    end
                    flag = false;
                    for i = 1:size(notPossible2,1)
                        if all(ismember(notPossible2(i,:),sta2))
                            flag = true;
                            break;
                        end   
                    end

                    if flag
                        continue
                    end

                    % write sub-configuration
                    subcon(nsub).nscan = 1;
                    subcon(nsub).scan(subcon(nsub).nscan).srcid = isrc1;
                    subcon(nsub).scan(subcon(nsub).nscan).startmjd = 0.0;
                    subcon(nsub).scan(subcon(nsub).nscan).nsta = pmcb2(ip).downstanum1;
                    for ista = 1 : subcon(nsub).scan(subcon(nsub).nscan).nsta
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).staid = pmcb2(ip).downstasn1(ista);
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).slewtime = 0.0;
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).startmjd = 0.0;
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).az       = 0.0;
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).el       = 0.0;
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).ha       = 0.0;
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).dc       = 0.0;
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).duration = 0.0;
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).endmjd   = 0.0;
                    end

                    subcon(nsub).nscan = 2;
                    subcon(nsub).scan(subcon(nsub).nscan).srcid = isrc2;
                    subcon(nsub).scan(subcon(nsub).nscan).startmjd = 0.0;
                    subcon(nsub).scan(subcon(nsub).nscan).nsta = pmcb2(ip).downstanum2;
                    for ista = 1 : subcon(nsub).scan(subcon(nsub).nscan).nsta
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).staid = pmcb2(ip).downstasn2(ista);
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).slewtime = 0.0;
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).startmjd = 0.0;
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).az       = 0.0;
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).el       = 0.0;
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).ha       = 0.0;
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).dc       = 0.0;
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).duration = 0.0;
                        subcon(nsub).scan(subcon(nsub).nscan).sta(ista).endmjd   = 0.0;
                    end
                    nsub = nsub + 1;
                end %for ip = 1 : np2
            end
        end
    end
end
