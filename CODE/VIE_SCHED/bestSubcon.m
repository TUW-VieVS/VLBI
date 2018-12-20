% Purpose  
%   Calculate the best subconfiguration(s)
%
%   INPUT:
%   source based:
%       bestSubcon( 'source',station,twin,source,srcat,catpair,staobs,srcobs,icatpair,srcn,obsmode,jpl,PARA );
%   station based:
%       bestSubcon( 'station',station,twin,source,staobs,srcobs,obsmode,PARA,starsource );
%
%   RETURN:
%   source based:  
%       returns best subcon
%   station based: 
%       returns array of possible subcons
%
%
% History  
%   2016-08-29   Matthias Schartner:  Created
%   2016-09-07   Matthias Schartner:  changes for singleCheckSource
%   2016-09-08   Matthias Schartner:  bugfix for source based scheduling
%   2016-11-03   Matthias Schartner:  singleCheckSource now saves az el ha dc in subcon
%   
function [ subcon, icatpair, srcn ] = bestSubcon( type,varargin )

%% source Based
if isequal(type,'source')
    station = varargin{1};
    twin = varargin{2};
    source = varargin{3};
    srcat = varargin{4};
    catpair = varargin{5};
    staobs = varargin{6};
    srcobs = varargin{7};
    icatpair = varargin{8};
    srcn = varargin{9};
    obsmode = varargin{10};
    jpl = varargin{11};
    PARA = varargin{12};
    
    loopnum = 0;
    srcnum = length(source);
    while true
        loopnum = loopnum + 1;
        if(loopnum > srcnum)
            subcon.nscan = 0;
            break;
        end

        [subcon1, icatpair, srcn] = ssubcon_src(station, twin, source, srcat, catpair, staobs, srcobs, icatpair, srcn, PARA); 

        if (subcon1.nscan == 0)
            subcon.nscan = 0;
            break;
        end
        % calculate start time
        [subcon2] = sstartmjd(source, station, twin, staobs, subcon1, PARA);
        if (subcon2.nscan == 0)
            continue;
        end
        % calculate duration
        [subcon] = sduration(source, station, twin, obsmode, subcon2, PARA);
        if (subcon.nscan  == 0)
            continue;
        end

        % check the network size
        netsize = 0;
        for iscan = 1 : subcon.nscan
            netsize = netsize + subcon.scan(iscan).nsta;
        end
        if (netsize < PARA.MIN_STANUM)
            continue;
        end
        
        % check for source structure study
        iforsi = 1;
        if (PARA.FORSI == 1)
            for iscan = 1 : subcon.nscan
                [obs_ok(iscan), grout(iscan)]= calc_soustruc_pscan(subcon.scan(iscan),source,station);
            end
            iforsi = sum(obs_ok(1:subcon.nscan));
        end
        if (iforsi == 0)
            continue;
        end

        [flag,~,subcon] = singleCheckSource(source, station, subcon, PARA, jpl);
        if flag == 1
            break
        end
    end
end
%% station Based
if isequal(type,'station')
    station = varargin{1};
    twin = varargin{2};
    source = varargin{3};
    staobs = varargin{4};
    srcobs = varargin{5};
    obsmode = varargin{6};
    PARA = varargin{7};
    starsource = varargin{8};
    
    subconall = ssubcon_sta(source, srcobs, station, staobs, PARA, starsource); 
    % check for source structure study
    if (PARA.FORSI == 1)
        for i = 1 : length(subconall)
            for iscan = 1 : subconall(i).nscan
                 subconall(i).scan(iscan).grout = 0;
            end
        end
    end

    % #### 1st sort ####
    % Sort all sub-cons calculated with "ssubcon_sta" initially => Then consider only the beste "PARA.SORTNUM" configurations
    [y] = ssubsort(source, station, staobs, srcobs, subconall, PARA); % y = sort index
    bestsubid = length(subconall) + 1; % index
    num = 0; % Counter for sub cons => Break condition: num == PARA.SORTNUM

    % #### Loop over sub-cons ####
    while true
        bestsubid = bestsubid - 1;
        if (bestsubid == 0) % Break condition: All sub-cons have been treated
            break;
        end

        % Apply sort index y and pick individual sub-cons, beginning with the best rated one.
        subcon1 = subconall(y(bestsubid)); 

        % ### calculate start time (1) ###
        [subcon2] = sstartmjd(source, station, twin, staobs, subcon1, PARA);
        if (subcon2.nscan < 1) % If all scans are thrown out => Discard the sub-con
            continue;
        end
        % ### calculate duration (1) ###
        [subcon3] = sduration(source, station, twin, obsmode, subcon2, PARA);
        if (subcon3.nscan  < 1) % If all scans are thrown out => Discard the sub-con
            continue;
        end
        subcon = subcon3;        

        % ### check the network size ###
        netsize = 0;
        for iscan = 1 : subcon.nscan
            netsize = netsize + subcon.scan(iscan).nsta;
        end
        if (netsize < PARA.MIN_STANUM)
            continue;
        end

        % ### check for source structure study ###
        iforsi = 1;
        if (PARA.FORSI == 1)
            for iscan = 1 : subcon.nscan
                [obs_ok(iscan), grout(iscan)]= calc_soustruc_pscan(subcon.scan(iscan),source,station);
                subcon.scan(iscan).grout = grout(iscan);
            end
            iforsi = sum(obs_ok(1:subcon.nscan));
        end
        if (iforsi == 0)
            continue;
        end

        % ### Save sub-con ###
        num = num + 1; % Increase the sub-con counter => break condition!
        subconall1(num) = subcon; % Save the individual sub-con which was currently treated to an array
        clear subcon;

        % ### Break condition: max. number of considered sub-cons) ###
        if (num == PARA.SORTNUM)   % PARA.SORTNUM = Num of subconfigurations found with station-based strategy for further consideration...
            break;
        end
    end % while true
    
    clear subcon;
    
    if num==0
        subcon.nscan = 0;
    else
        subcon = subconall1;
    end
    
    % icatpair and srcn not needed for station based
    icatpair = []; 
    srcn = [];
end

end

