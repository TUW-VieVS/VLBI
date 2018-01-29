% Purpose  
%   Sort all the possible subnets by the source structure study, num of obs, sky coverage, end time, source weight.
% History  
%   2010-04-06   Jing SUN   Created
%
% Changes
%   2015-07-15, A. Hellerschmied: Comments added.
%   2015-11-02, A. Hellerschmied: Comments added.
%   2016-06-17, M. Schartner: changes to improve Speed, 
%                             changes for vectorised zazel_s.m function
%   2016-09-06, M. Schartner: "score" as 2nd return value

function [yy, xysort] = ssubsort(source, station, staobs, srcobs, subcon, PARA)

% ##### Loop over all possible sub-cons #####
% Calculate all values which are used as "wheight factors" for the sorting
nmax = length(subcon);
subss = zeros(nmax,1);
substasky = zeros(nmax,1);      % sky coverage
subobsnum = zeros(nmax,1);      % 
subendtim = zeros(nmax,1);
subsweight = zeros(nmax,1);      % Jing Sun, 2014 Aug 13 - source weight
subsrcscan = zeros(nmax,1);      % Jing Sun, 2014 Aug 13

llh = [station.llh];
lon = llh(1:3:end);
lat = llh(2:3:end);

for i = 1 : nmax 
    % #### Loop over all scans (1 or 2) within one sub-con ####
    for j = 1 : subcon(i).nscan 
        
        if (PARA.FORSI == 1) % Source structure
            subss(i) = subss(i) + subcon(i).scan(j).grout;
        end
        
        srcid = subcon(i).scan(j).srcid;
        subsweight(i) = subsweight(i) + source(srcid).weight;     % Jing Sun, 2014 Aug 13
        subsrcscan(i) = subsrcscan(i) + srcobs(srcid).nscan;      % Jing Sun, 2014 Aug 13
        
        
        staid_all = [subcon(i).scan(j).sta.staid];
        [az, el] = zazel_s([staobs(staid_all).endmjd], lon(staid_all), lat(staid_all), source(srcid).ra, source(srcid).de); % Az/El calculation
        
        % #### Loop over all stations in this scan ####
        for k = 1 : subcon(i).scan(j).nsta  
            
            % ### Calc. sky coverage number ###
            staid = subcon(i).scan(j).sta(k).staid;
            [cat] = zcoverage(az(k), el(k)); % Calculate sky coverage number at the station.  
            ifnewcat = 1;  % New scan for sky coverage number? yes ==> =1
            
            % ### Check if the sky coverage number is new (the sky coverage increases) ###
            for icat = 1 : staobs(staid).ncat
                if(staobs(staid).catsn(icat)==cat)
                    ifnewcat = 0;
                    break;
                end
            end
            if (ifnewcat == 1)
                substasky(i) = substasky(i) + 1; % If new number => Increase sky coverage count (max. = number of stations station)
            end
            
            % ### Get the latest "observation end time" of all stations within all scans of the sub-con ###
            % ==> end time of the whole sub-con
            if (subcon(i).scan(j).sta(k).endmjd > subendtim(i)) 
                subendtim(i) = subcon(i).scan(j).sta(k).endmjd;
            end
        end
        
        % ### Number of observations ###
        subobsnum(i) = subobsnum(i) + subcon(i).scan(j).nsta * (subcon(i).scan(j).nsta - 1)/2;
    end 
    
    % ### source weight ###
    subsweight(i) = subsweight(i) / subcon(i).nscan; % David Mayer, 2014 Oct 28 - norm factor -> otherwise a subcon with more subnets is weighted better
    
end % for i = 1 : length(subcon) 




% #################################################################
% de/increase the opportunity for scheduling
% PARA.UPSTA, PARA.DOWNSTA
% Add if-routine here!!!!
stanum = length(station);
upstaid = 0;
downstaid = 0;
for ista = 1 : stanum
    if strcmp(deblank(PARA.UPSTA), deblank(station(ista).name))
        upstaid = ista;
        break;
    end
end
for ista = 1 : stanum
    if strcmp(deblank(PARA.DOWNSTA), deblank(station(ista).name))
        downstaid = ista;
        break;
    end
end

for i = 1 : length(subcon)
    for j = 1 : subcon(i).nscan
        for k = 1 : subcon(i).scan(j).nsta
            staid = subcon(i).scan(j).sta(k).staid;
            if (staid == upstaid)
                substasky(i) = substasky(i) + 1; % Increase the sky coverage number
            elseif (staid == downstaid)
                if (substasky(i) >= 1)
                    substasky(i) = substasky(i) - 1; % Decrease the sky coverage number
                else 
                    substasky(i) = 0;
                end
            end
        end
    end
end
% #################################################################

% ##### Calc weights for all sub-cons ##### 
wnos =  PARA.WEIGHT_NUMBER_OF_OBS;
wsky =  PARA.WEIGHT_SKY_COVERAGE;
wset =  PARA.WEIGHT_SCAN_END_TIME;

for i = 1 : length(subcon) % all sub cons
    
    % Init.:
    xysort(i) = 0;
    
    if (PARA.FORSI == 1)                                        % Source Structure
        xysort(i) = xysort(i) + subss(i)*10;
    end
    
    xysort(i) = xysort(i) + subobsnum(i)/(stanum*(stanum-1)/2)*wnos; % Number of observations (range of values: 0 - 1)
    xysort(i) = xysort(i) + substasky(i)/stanum*wsky;                % Sky coverage (sum of all stations) (range of values: 0 - 1)
    
    maxendtim = max(subendtim(1:length(subcon)));               
    minendtim = min(subendtim(1:length(subcon)));
    if ((maxendtim-minendtim) > 0.0d0)                          % Scan end time
        xysort(i) = xysort(i) + ((maxendtim-subendtim(i))/(maxendtim-minendtim))*wset;   
    end
    
    xysort(i) = xysort(i) * subsweight(i);      % David Mayer, 2014 Oct 28 , source weight
    
    if ((max(subsrcscan(:))-min(subsrcscan(:))) > 0.0d0) && PARA.DISTRIBUTE % David Mayer, 2014 Sep 3 
        xysort(i) = xysort(i) + subsweight(i);      % Jing Sun, 2014 Aug 13 
        xysort(i) = xysort(i) + (max(subsrcscan(:))-subsrcscan(i))/(max(subsrcscan(:))-min(subsrcscan(:)));      % Jing Sun, 2014 Aug 13
    end
end

% ##### Sort #####
[xx, yy] = sort(xysort); % yy = sort indices; xx = sorted xysort vector


