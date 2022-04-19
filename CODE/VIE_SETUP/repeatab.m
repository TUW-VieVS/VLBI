% Calculate baseline length repeatability for longer time-spans, ie
% including breaks
%
% ==================================================================
% Main VieVS function for calculating baseline length repeatability! 
% Should provide all flexibility
%  - simple plot
%  - text output
%  - variable output
% ==================================================================
%
% INPUT (all are varargin)
%  subf             subfolder of LEVEL3/ and LEVEL1/ folders. 
%                   Default: ''
%  process_list     process_list (either variable or mat file)
%                   Default: 'process_list.mat' in WORK directory
%  limitation       minimum number of observations of baselines per
%                   "break-time-span".
%                   Default: 10
%  outfile          full out path (incl fname). Is taken as is: If not path
%                   is given, file is put to Matlab's current folder.
%                   Default: [] -> no outfile
%  basOutFname      filename for bas out (might be needed in some cases)
%  makeFig          1 = plot. 
%                   Default: 0
%  printToCommand   1 = print baselines to command window. 
%                   Default: 0
%  superstations    superstations file (either variable or fname). 
%                   Default: loaded from TRF folder (only works if current
%                   folder is WORK folder)
%  rigFormErr       Rigorous treatment of formal errors. 1=rigorous,
%                   0=simple calculation. Default: 1;
%  allFiles         cell, containing 5 cells of loaded files
%                   {1} x_files    struct (length = nSessions) containing x_ files. 
%                       Needed when files are already loaded (eg VieVS)
%                       [] -> not given
%                   {2} ant_files  struct (length = nSessions) containing antenna files
%                   {3} atpa_files struct (length = nSessions) containing atpa_ files
%                   {4} opt_files  struct (length = nSessions) containing apt_ files
%
% OUTPUT (all are varargout)
%  blr              baseline length repeatability (m) with limit=10 (min
%                   numObs of baseline per "break-time-span". For more
%                   flexibility, use following output parameters! Baselines
%                   less than limitations are NaN.
%  wblr  			weighted blr (according formal errors)
%  bl               baseline lengths (m)
%  blnames          baseline names
%  blr_all          baseline length repeatabilites (m) for all
%                   break-time-spans. Matrix (numBaselines x
%                   maxNumOfBreak-Time-Spans). Other entries (than needed)
%                   are NaN!
%  wblr_all         weighted blr_all
%  nobs_all         Number of observation. Corresponding to blr_all. 
%                   Other fields are NaN.
%  mjds_all         mean mjds of observations. Corresponding to 
%                   blr_all. Other fiels are NaN.
%  mjdStart_all     Min mjd of observations. Corresponding to 
%                   blr_all. Other fiels are NaN.
%  mjdEnd_all       Max mjd of observations. Corresponding to 
%                   blr_all. Other fiels are NaN.
%
% LOG
%  30.11.2015 M. Madzak: Created
%
% CHANGES
%  2016-01-11 M. Madzak: Bug fix
%  2016-02-16 A. Hofmeister: Bug fix
%  2016-02-26 M. Madzak: Added baseline length estimate even if no repeatab can be calculated (only one estimate of that baseline)
%  2016-03-22 A. Hofmeister: change degree of freedom for wblr calculation in order to have an
%                            unbiased wblr,
%                            use inverse of squared formal baseline length errror as weight,
%                            change determination of blr, wblr, mjd and bl from normal mean to a
%                            weighted mean estimation of the individual per break values
% 2016-09-21, A. Girdiuk:  warning messages are added in case of absence of station coordinates estimates
% 2017-05-31, A. Hellerschmied: costants.m added.
%

function varargout = repeatab(varargin)


%% default input arguments
subf='';
process_list='process_list.mat';
limitation=10;      % minimum number of baseline observations
outfile=[];
basOutFname=[];
makeFig=0;
printToCommand=0;
superstations='../TRF/superstation.mat';
rigFormErr=1;

% Add constants (global variables)
% Note: constants.m is located in VIE_MOD_xx
constants;


%% get input arguments
if nargin>0
    if ~isempty(varargin{1}); subf=varargin{1}; end
    if nargin>1
        process_list=varargin{2};
        if nargin>2
            if ~isempty(varargin{3}); limitation=varargin{3}; end
            if nargin>3
                outfile=varargin{4};
                if nargin>4
                    basOutFname=varargin{5};
                
                    if nargin>5
                        if ~isempty(varargin{6}); makeFig=varargin{6}; end
                        if nargin>6
                            if ~isempty(varargin{7}); printToCommand=varargin{7}; end
                            if nargin>7
                                if ~isempty(varargin{8}); superstations=varargin{8}; end
                                if nargin>8
                                    if ~isempty(varargin{9}); rigFormErr=varargin{9}; end
                                end
                            end

                        end
                    end
                
                end
            end
        end
    end
end
 
% get files (x_,...)
x_filesGiven=0;
x_files=[];
antFilesGiven=0;
antFiles=[];
atpaFilesGiven=0;
atpaFiles=[];
optFilesGiven=0;
optFiles=[];
if nargin>9
    if ~isempty(varargin{10})
        if ~isempty(varargin{10}{1})
            x_files=varargin{10}{1};
            x_filesGiven=1;
        end
        if ~isempty(varargin{10}{2})
            antFiles=varargin{10}{2};
            antFilesGiven=1;
        end
        if ~isempty(varargin{10}{3})
            atpaFiles=varargin{10}{3};
            atpaFilesGiven=1;
        end
        if ~isempty(varargin{10}{4})
            optFiles=varargin{10}{4};
            optFilesGiven=1;
        end
    end        
end

% % load it if it's a filename
% if strcmpi(process_list(end-3:end),'.mat') % if it's a filename
%     process_list=load(process_list); 
%     fpl=fieldnames(process_list); process_list=process_list.(fpl{1});
% end

%% Get baselines
bas=bas_out(process_list,subf,basOutFname,...
    x_files,antFiles,atpaFiles,optFiles,rigFormErr);

if isempty([bas.ap])
    varargout{1}=nan;
    varargout{2}=nan;
    varargout{3}=nan;        % baseline lengths (m)
    varargout{4}={};% baseline names
    fprintf('I can not calculate baselines repeatability because bas_out.m output was empty\n');
    return
end

% load superstations file
if strcmpi(superstations(end-3:end),'.mat') % if it's a filename
    superstations=load(superstations);
    fieldsSuper=fieldnames(superstations); 
    superstations=superstations.(fieldsSuper{1});
end


%% preallocate
bl=ones(length(bas),1)*NaN;
mjd=bl;
blr=bl;
wblr=bl; %weighted RMS

bl_all=ones(length(bas),99)*NaN; % full out matrix -> row=baseline, col=one "break-timespan"
blr_all=bl_all;
wblr_all=bl_all;
nobs_all=bl_all;
mjds_all=bl_all; % mean of mjds
mjdStart_all=bl_all;
mjdEnd_all=bl_all;


%% loop over baselines
for k=1:length(bas)

    % get station indices in superstation file
    index_1=find(strcmpi({superstations.name}, bas(k).name(1:8)));
    index_2=find(strcmpi({superstations.name}, bas(k).name(10:17)));
    
    if isempty(index_1)
        curS_trim=strtrim(bas(k).name(1:8));
        curS__=strrep(curS_trim,' ', '_'); % might be only 6 chars
        curS=[curS__, repmat(' ', 1, 8-length(curS__))];
        index_1=find(strcmpi({superstations.name}, curS));
    end
    if isempty(index_2)
        curS_trim=strtrim(bas(k).name(10:17));
        curS__=strrep(curS_trim,' ', '_'); % might be only 6 chars
        curS=[curS__, repmat(' ', 1, 8-length(curS__))];
        index_2=find(strcmpi({superstations.name}, curS));
    end
     
    % get start/end epochs of breaks for both indices
    if ~isfield(superstations(index_1).vievsTrf.break,'start') % some stations have no start field
        i1_breakStart=0;
        i1_breakEnd=999999;
    else
        i1_breakStart=[superstations(index_1).vievsTrf.break.start];
        i1_breakEnd=[superstations(index_1).vievsTrf.break.end];
    end
    if ~isfield(superstations(index_2).vievsTrf.break,'start') % some stations have no start field
        i2_breakStart=0;
        i2_breakEnd=999999;
    else
        i2_breakStart=[superstations(index_2).vievsTrf.break.start];
        i2_breakEnd=[superstations(index_2).vievsTrf.break.end];
    end
    
    % get break index for all bas-estimates
    allObsEp=bas(k).mjd; % for all those epochs -> get breaks of both (participating stations)
    breakInd=zeros(2,length(allObsEp)); % contains breaks for both stations (row) for all epochs (columns) of bas-observations
 
    for iEp=1:length(allObsEp)
        
        breakInd(1,iEp)=find(i1_breakStart<=allObsEp(iEp) & ...
            i1_breakEnd>allObsEp(iEp));
        
        breakInd(2,iEp)=find(i2_breakStart<=allObsEp(iEp) & ...
            i2_breakEnd>allObsEp(iEp));
    end
    allBreakCombis=unique(breakInd','rows')';
    nUniqBreaks=size(allBreakCombis,2);
    blrPerBreak=ones(1,nUniqBreaks)*NaN;
    wblrPerBreak=blrPerBreak;
    mjdPerBreak= blrPerBreak;
    meanPerBreak=blrPerBreak;
    for kUniqBreak=1:nUniqBreaks
        % get observations of current unique combination
        curObs=breakInd(1,:)==allBreakCombis(1,kUniqBreak) & ...
            breakInd(2,:)==allBreakCombis(2,kUniqBreak);
        curBL=mean(bas(k).corr(curObs));
        curMjd=mean(bas(k).mjd(curObs));
        mjdStart_all(k,kUniqBreak)=min(bas(k).mjd(curObs));
        mjdEnd_all(k,kUniqBreak)=max(bas(k).mjd(curObs));
        curBLR=NaN; % intialization is necessary for each new baseline cycle
        curWBLR=NaN; % intialization is necessary for each new baseline cycle
        
        if sum(curObs)>1 % to distinguish equal lengths (->std=0) from one observation (-> std remains NaN since curBLR and curWBLR initialized as NaN)
            
            curMjds=bas(k).mjd(curObs); % sort by mjd...
            curCorr=bas(k).corr(curObs);
			curWeights=1./(bas(k).mbas(curObs).^2); % can be Inf (if formal error is 0), inverse of squared formal error
            % treat infinite numbers
            if sum(isinf(curWeights))>0
                curWeights(~isinf(curWeights))=0;
                curWeights(isinf(curWeights))=1;
            end
            [curMjds_sort,sI]=sort(curMjds);
            curCorr_sort=curCorr(sI);
			curW_sort=curWeights(sI);
            p_linFit=polyfit(curMjds_sort,curCorr_sort,1);
            
            
            curCorr_sort_linFitRem=curCorr_sort-polyval(p_linFit,curMjds_sort); 
            curBLR=std( curCorr_sort_linFitRem );
            
            sum_curW_sort=sum(curW_sort); % sum of weights
			curWBLR=sqrt( sum(curW_sort.*(curCorr_sort_linFitRem-mean(curCorr_sort_linFitRem)).^2) / (sum_curW_sort - sum(curW_sort.^2)/sum_curW_sort) ); % unbiased wblr estimation
            
        end        
        
        if sum(curObs)>=limitation % if there are enough observations inside current break

            blrPerBreak(1,kUniqBreak)=curBLR;
			wblrPerBreak(1,kUniqBreak)=curWBLR;
            mjdPerBreak(1,kUniqBreak)=curMjd;

        end
		meanPerBreak(1,kUniqBreak)=curBL; % save mean BL even if there is just 1 (or <limitation) estimates of that baseline
		
        % save for later use
        blr_all(k,kUniqBreak)=curBLR;
		wblr_all(k,kUniqBreak)=curWBLR;
        bl_all(k,kUniqBreak)=curBL;
        nobs_all(k,kUniqBreak)=sum(curObs);
        mjds_all(k,kUniqBreak)=mean(bas(k).mjd(curObs));
    end

    % Get the final blr, wblr, mjd for one baseline: weighted mean over all proper (more than limit
    % observations) breaks.
    % A valid bl is determined also if limit is not reached by the number of observations.
    % The values are determined as weighted mean from the values of all individual break
    % intervals. The weights are the number of observations in each break interval.
    % In case not any break interval delivered a valid value, if the number of observations did not
    % reach the limit, the result is set to NaN.
    
    isvalid_blrPerBreak= ~isnan(blrPerBreak);
    if any(isvalid_blrPerBreak)
        blr(k,1)= sum( nobs_all(k,isvalid_blrPerBreak) .* blrPerBreak(isvalid_blrPerBreak) ) / sum( nobs_all(k,isvalid_blrPerBreak) );
    else
        blr(k,1)= NaN;
    end
    
    isvalid_wblrPerBreak= ~isnan(wblrPerBreak);
    if any(isvalid_wblrPerBreak)
        wblr(k,1)= sum( nobs_all(k,isvalid_wblrPerBreak) .* wblrPerBreak(isvalid_wblrPerBreak) ) / sum( nobs_all(k,isvalid_wblrPerBreak) );
    else
        wblr(k,1)= NaN;
    end
    
    isvalid_mjdPerBreak= ~isnan(mjdPerBreak);
    if any(isvalid_mjdPerBreak)
        mjd(k,1)= sum( nobs_all(k,isvalid_mjdPerBreak) .* mjdPerBreak(isvalid_mjdPerBreak) ) / sum( nobs_all(k,isvalid_mjdPerBreak) );
    else
        mjd(k,1)= NaN;
    end
    
    isvalid_meanPerBreak= ~isnan(meanPerBreak);
    if any(isvalid_meanPerBreak)
        bl(k,1)= sum( nobs_all(k,isvalid_meanPerBreak) .* meanPerBreak(isvalid_meanPerBreak) ) / sum( nobs_all(k,isvalid_meanPerBreak) );
    else
        bl(k,1)= NaN;
    end
    
end

% delete not needed cols (due to preallocation)
cols2del=sum(isnan(nobs_all))==size(nobs_all,1);
blr_all(:,cols2del)=[];
wblr_all(:,cols2del)=[];
%bl_all(:,cols2del)=[]; % not output, thus commented
nobs_all(:,cols2del)=[];
mjds_all(:,cols2del)=[];
mjdStart_all(:,cols2del)=[];
mjdEnd_all(:,cols2del)=[];

%% OUTPUT

% 1. output to textfile
if ~isempty(outfile)
    curClock=clock;
    fid=fopen(outfile, 'w');
    fprintf(fid,['# Baseline length repeatability from Vienna VLBI Software\n',...
        '# Mean values over all break-time-spans\n',...
        '# Linear fit removed before calculation of standard deviation\n#\n',...
        '# Breaks (Earthquakes) are taken from the vievsTrf\n#\n',...
        '# Created on %02.0f.%02.0f.%04.0f %02.0f:%02.0f:%02.0f\n',...
        '# by function repeatab.m\n#\n',...
        '# Min number of observations of baseline (in every break-time-span): %1.0f\n',...
        '# Number of baselines: %1.0f\n#\n',...
        '# col1 (cols 01-17)  baseline name \n',...
        '# col2 (cols 19-28)  mean epoch (mjd)\n',...
        '# col3 (cols 30-42)  mean baseline length in meters\n',...
		'# col4 (cols 44-49)  baseline length repeatability in cm\n',...
        '# col5 (cols 51-56)  weighted baseline length repeatability in cm\n#\n',...
        '# col6 (cols 58-62)  number of sessions\n#\n'],...
        curClock(3:-1:1), curClock(4:6), limitation,...
        sum(~isnan(blr) & blr~=0));
    
    for iB=1:length(bas)
        if ~isnan(blr(iB)) && blr(iB)~=0 % if there was more than one estimate
            fprintf(fid, '%-17s %10.4f %13.4f %6.2f %6.2f %5.0f\n',...
                bas(iB).name, mjd(iB), bl(iB), blr(iB)*100, wblr(iB)*100, nobs_all(iB));
        end
    end    
    fclose(fid);    
end

% 2. create figure
if makeFig==1
    figure('name', 'Baseline length repeatability');
    plot(bl/1000,blr*100,'ko', 'markerfacecolor', 'k');
    title('(Unweighted) baseline length repeatability');   
    xlabel('Baseline length (km)');
    ylabel('Standard deviation (cm)');
end

% 3. print to command windows
if printToCommand==1
    fprintf(' nr baselineName----- basLength--- BLR---- weighted\n');
    % print all found baselines (even those with no observatios... they
    % will print NaN!)
    for iB=1:length(bas)
        fprintf('%3.0f %s %11.2fm %5.2fcm %5.2fcm\n', ...
            iB, bas(iB).name, bl(iB), blr(iB)*100, wblr(iB)*100);
    end    
end


%% Assign variable outputs
varargout{1}=blr;
varargout{2}=wblr;
varargout{3}=bl;        % baseline lengths (m)
varargout{4}={bas.name}';% baseline names
varargout{5}=blr_all;   % baseline length repeatabilites (m)
varargout{6}=wblr_all;   % weighted baseline length repeatabilites (m)
varargout{7}=nobs_all;  % Number of observation.
varargout{8}=mjds_all;  % mean mjds
varargout{9}=mjdStart_all;
varargout{10}=mjdEnd_all;

