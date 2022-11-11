% This function reads sinex input
% First version: Only SITE/ID, SOLUTION/ESTIMATE blocks
%
%
% INPUT
%  varargin -> each 'parameter' consists of two parts:
%              general: 'PropertyName',PropertyValue
%              - 'blocks2read'    'all' or specify block(s) to be read in
%                                 (char array or cell array)
%
% OUTPUT (varargout)
%  blockList   cellstr of blocks (if 'blocks2read' given, this list is
%              equal to that input
%  outBlocks   cell (one entry = one block) containing the data of the
%              blocks. If 'all' blocks are read, the order is as in the snx
%              file. If 'blocks2read' was specified, the order is according
%              to that list. In any case, the order is given in
%              'blockList'. (Blocks from 'blocks2read' not existing in snx 
%              file will have an empty entry in cell 'outBlocks')
%
% CREATED: Matthias Madzak, 2016-05-10
%
% some infos
%  the site-codes are corresponding in SITE/ID, SOLUTION/EPOCHS and
%  SOLUTION/ESTIMATES ... -> same indices mean same site-code!

function [varargout] = readSnx(infile, varargin)

% infile='../TRF/data/DTRF2020P_VLBI.snx';
% 
% varargin{1}='blocks2read';
% varargin{2}={'SOLUTION/EPOCHS', 'SITE/ID'};

%% GET OPTIONS
% get defaults
readAllBlocks=1;
readBlocks=''; 

if nargin>1
    curVarArgNr=1;
    % varargin must be even number
    if mod(length(varargin),2)~=0
        fprintf('input to readSnx must contain ''propertyNames'' and property value (ie even number of options!)\n');
        return
    end
    
    % as long as there are more options in varargin
    while length(varargin)>curVarArgNr
        switch varargin{curVarArgNr}
            case 'blocks2read'
                curVarArgNr=curVarArgNr+1;
                if strcmpi(varargin{curVarArgNr},'all')
                    readAllBlocks=1;
                else
                    readAllBlocks=0;
                    readBlocks=cellstr(varargin{curVarArgNr});
                    nReadBlocks=length(readBlocks);
                end
            otherwise
                fprintf('''PropertyName'' %s is unknown! function readSnx.m\n',...
                    varargin{curVarArgNr});
        end
    end

    
end

%% "GLOBAL" preallocation
nPreallSites=5000;
siteCodesGlobal=zeros(nPreallSites,1); % do this "globally" as site/id block need not be read in before other blocks also requiring the codes

%% READ FILE


fid=fopen(infile);

blocksReadDoneInd=0;
while ~feof(fid)
    curl=fgetl(fid);
    if isempty(curl)
        curl=' ';
    end
    
    % if new block
    if strcmpi(curl(1),'+')
        curBlockName=deblank(curl(2:end));
        if readAllBlocks || sum(strcmpi(curBlockName,readBlocks))>0
            
            curBlockDone=0;
            if strcmp(curBlockName,'SITE/ID')
                outBlock{blocksReadDoneInd+1,1}=readSiteIdBlock(fid);
                curBlockDone=1;
            elseif strcmp(curBlockName,'SOLUTION/EPOCHS')
                outBlock{blocksReadDoneInd+1,1}=readSolutionEpochsBlock(fid);
                curBlockDone=1;
            elseif strcmp(curBlockName,'SOLUTION/ESTIMATE')
                outBlock{blocksReadDoneInd+1,1}=readSolutionEstimateBlock(fid);
                curBlockDone=1;
            else
                fprintf('block ''%s'' not done yet!\n',...
                    curBlockName);
            end
            
            if curBlockDone
                blocksReadNames{blocksReadDoneInd+1}=curBlockName;
                blocksReadDoneInd=blocksReadDoneInd+1;
            end
        end
    end
end

fclose(fid);



% "SORT" OUTPUT (if block names were specified in input)
if readAllBlocks==0
    takeBlocksInd=zeros(length(readBlocks),1);
    outBlockFinal=cell(length(readBlocks),1);
    blockList=cell(length(readBlocks),1); % names of blocks for final output
    for kBlockGiven=1:length(readBlocks)
        curLog=strcmpi(readBlocks{kBlockGiven},blocksReadNames);
        if sum(curLog)>0 % if the block was read in (ie found in snx file)
            takeBlocksInd(kBlockGiven,1)=find(curLog);
        end
    end
    outBlockFinal(takeBlocksInd~=0,1)=outBlock(takeBlocksInd(takeBlocksInd~=0));
    blockList(takeBlocksInd~=0,1)=blocksReadNames(takeBlocksInd(takeBlocksInd~=0));
else
    outBlockFinal=outBlock;
    blockList=blocksReadNames;

end

% assigne output
if nargout>0
    varargout{1}=blockList;
    varargout{2}=outBlockFinal;
end

aaa=5;


% --------------------------------------------------------------

function solEp=readSolutionEpochsBlock(fid)
    deleteEmptyBreaks=1; % "wrong" sinex (solutino number 1, 5 exist, others not)
    curBlockName='SOLUTION/EPOCHS';
    curl=fgetl(fid); if isempty(curl); curl=' '; end
    
    solEp=struct('code', [], 'pt', [], 'break', []);
    while ~strcmpi(deblank(curl),['-',curBlockName])
        % if not a comment line
        if ~strcmpi(curl(1),'*')
            solNr=str2double(curl(10:13)); 
            % if we already have that code (in the "global" vector)
            if sum(siteCodesGlobal==str2double(curl(2:5)))>0
                curCodeInd=find( siteCodesGlobal==str2double(curl(2:5)) );
                
%                 if length(solEp)>=curCodeInd
%                     if ~isempty(solEp(curCodeInd).break)
%                         solNr=length(solEp(curCodeInd).break)+1;
%                     end
%                 end % Sometimes the sol-nr in the sinx file is wrong/bad!
                
            else
                curCodeInd=length(siteCodesGlobal~=0)+1;
                siteCodesGlobal(curCodeInd)=str2double(curl(2:5));
            end
            
            solEp(curCodeInd).code=str2double(curl(2:5));
            solEp(curCodeInd).pt=curl(7:8);
            solEp(curCodeInd).t=curl(15);
            
            curYr=str2double(curl(17:18));
            curDoy=str2double(curl(20:22));
            curSecod=str2double(curl(24:28));
            solEp(curCodeInd).break(solNr).start=date2mjd([curYr 1 1])+...
                curDoy-1+curSecod/60/60/24;
            
            curYr=str2double(curl(30:31));
            curDoy=str2double(curl(33:35));
            curSecod=str2double(curl(37:41));
            solEp(curCodeInd).break(solNr).end=date2mjd([curYr 1 1])+...
                curDoy-1+curSecod/60/60/24;
            
            curYr=str2double(curl(43:44));
            curDoy=str2double(curl(46:48));
            curSecod=str2double(curl(50:54));
            solEp(curCodeInd).break(solNr).meanEpoch=date2mjd([curYr 1 1])+...
                curDoy-1+curSecod/60/60/24;
            
            % try to find current code in epoch struct
%             [ep.code]==inf;
%             code(kSite)=str2double(curl(2:5));
%             pt{kSite}=curl(7:8);
%             domes{kSite}=curl(10:18);
%             t{kSite}=curl(20);
%             station{kSite}=curl(22:29);
%             descr{kSite}=curl(31:43);
            % TO DO: Approx coordinates!
            
        end
        
        % read next line
        curl=fgetl(fid);
        if isempty(curl); curl=' '; end
        
    end
    
    % delete not wanted lines
%     if preallLength>kSite
%         code(kSite+1:end)=[];
%         pt(kSite+1:end)=[];
%         domes(kSite+1:end)=[];
%         t(kSite+1:end)=[];
%         station(kSite+1:end)=[];
%         descr(kSite+1:end)=[];
%     end

    % delete empty breaks (which are there if sinex is badly written)
    if deleteEmptyBreaks
        % for all available stations/epo
        for k=1:length(solEp)
            if ~isempty(solEp(k).break) %  this might be empty (even though stupid)
                allFields=fieldnames(solEp(k).break); % fields of break

                nBr=length(solEp(k).break);
                breaksEmpty=zeros(nBr,1);
                for kBr=1:nBr
                    nFieldsEmpty=0;
                    for kF=1:length(allFields)
                        if isempty(solEp(k).break(kBr).(allFields{kF}))
                            nFieldsEmpty=nFieldsEmpty+1;
                        end
                    end
                    if nFieldsEmpty==length(allFields)
                        breaksEmpty(kBr,1)=1;
                    end
                end
                % delete the empty breaks
                solEp(k).break(logical(breaksEmpty))=[];
            end
        end
    end
end

% --------------------------------------------------------------
function solEst=readSolutionEstimateBlock(fid)
    deleteEmptyBreaks=1; % "wrong" sinex (solutino number 1, 5 exist, others not)
    
    curBlockName='SOLUTION/ESTIMATE';
    curl=fgetl(fid); if isempty(curl); curl=' '; end
    
    statStruct=struct('code', [], 'break', [], 'pt', []);
    solEst=struct('stat', statStruct); % there could also be EOP or others
    
    allFieldsOfBreak={'x','x_sigma','y','y_sigma','z','z_sigma',...
        'vx','vx_sigma','vy','vy_sigma','vz','vz_sigma'};
    
    while ~strcmpi(deblank(curl),['-',curBlockName])
        % if not a comment line
        if ~strcmpi(curl(1),'*')
            curCodeStr=curl(15:18);
            
            if strcmpi(curCodeStr,'----')
                fprintf('Not done yet!\nkeyboard\nfunction readSnx.m\n');
                keyboard;
            else % we have a station dependent parameter
                solNr=str2double(curl(23:26));
                curCode=str2double(curCodeStr);
                if curCode==1857
%                     keyboard;
                end
                % if we already have that code (in the "global" vector)
                if sum(siteCodesGlobal==curCode)>0
                    curCodeInd=find( siteCodesGlobal==curCode );
                    solNr=str2double(curl(23:26)); %  does not work if
%                     sinex are badly (wronly) written
                    
%                      The following (now commented) block is if the
%                      solution number in SINEX file is bad/wrong...
%                      still, it is a bad solution...
%                     % if there is not even an index of that station
%                     if length(solEst.stat)<curCodeInd
%                         solNr=1;
%                     else
%                         if ~isfield(solEst.stat(curCodeInd),'break')
%                             solNr=1;
%                         else
%                             % if the break-field is empty
%                             if isempty(solEst.stat(curCodeInd).break)
%                                 solNr=1;
%                             else
%                                 solNr=length(solEst.stat(curCodeInd).break);
% 
%                                 allFieldsExist=1; % this remains one if alls fields exist and being not empty
%                                 for kField=1:length(allFieldsOfBreak)
%                                     % if the field(s) exist in current
%                                     % break
%                                     if ~isfield(solEst.stat(curCodeInd).break(solNr), allFieldsOfBreak{kField})
%                                         allFieldsExist=0;
%                                         break;
%                                     else
%                                         if isempty(solEst.stat(curCodeInd).break(solNr).(allFieldsOfBreak{kField}))
%                                             allFieldsExist=0;
%                                             break;
%                                         end
%                                     end
%                                 end
% 
%                                 if allFieldsExist==1
%                                     solNr=solNr+1;
%                                 end
%                             end
%                         end
%                     end
                        
             
                    
                else
                    curCodeInd=length(siteCodesGlobal~=0)+1;
                    siteCodesGlobal(curCodeInd)=curCode;
                    
                end
                
                % if we may access the current index
                if length(solEst.stat)>=curCodeInd
                    if isempty(solEst.stat(curCodeInd).code)
                        solEst.stat(curCodeInd).code=curCode;
                        solEst.stat(curCodeInd).pt=curl(20:21);
                        solEst.stat(curCodeInd).s=curl(46);
                    end
                else
                    solEst.stat(curCodeInd).code=curCode;
                    solEst.stat(curCodeInd).pt=curl(20:21);
                    solEst.stat(curCodeInd).s=curl(46);
                end
                
                curVal=str2double(curl(48:68));
                curStd=str2double(curl(70:80));
                
                % check unit (make good unit)
                curUnit=strtrim(curl(41:44));
                if ~(strcmpi(curUnit,'m') || strcmpi(curUnit,'m/y'))
                    fprintf('WARNING: Unit %s in snx file unknown - be sure to read in proper!!!\n',...
                        curUnit);
                end
                
                switch curl(8:11)
                    case 'STAX'
                        solEst.stat(curCodeInd).break(solNr).x=curVal;
                        solEst.stat(curCodeInd).break(solNr).x_sigma=curStd;
                        epochGiven=1;
                    case 'STAY'
                        solEst.stat(curCodeInd).break(solNr).y=curVal;
                        solEst.stat(curCodeInd).break(solNr).y_sigma=curStd;
                        epochGiven=1;
                    case 'STAZ'
                        solEst.stat(curCodeInd).break(solNr).z=curVal;
                        solEst.stat(curCodeInd).break(solNr).z_sigma=curStd;
                        epochGiven=1;
                    case 'VELX'
                        solEst.stat(curCodeInd).break(solNr).vx=curVal;
                        solEst.stat(curCodeInd).break(solNr).vx_sigma=curStd;
                        epochGiven=1;
                    case 'VELY'
                        solEst.stat(curCodeInd).break(solNr).vy=curVal;
                        solEst.stat(curCodeInd).break(solNr).vy_sigma=curStd;
                        epochGiven=1;
                    case 'VELZ'
                        solEst.stat(curCodeInd).break(solNr).vz=curVal;
                        solEst.stat(curCodeInd).break(solNr).vz_sigma=curStd;
                        epochGiven=1;
                    otherwise
                        fprintf('Type %s in SOLUTION/ESTIMATE block not defined yet \n',...
                            curl(8:11));
                        epochGiven=0;
                end
                        
                if epochGiven
                    curYr=yy2yyyy(str2double(curl(28:29)));
                    curDoy=str2double(curl(31:33));
                    curSecod=str2double(curl(35:39));
                    
                    solEst.stat(curCodeInd).break(solNr).epoch=...
                        date2mjd([curYr 1 1])+...
                        curDoy-1+curSecod/60/60/24;
                end
            end
            
        end
        % read next line
        curl=fgetl(fid);
        if isempty(curl); curl=' '; end
    end
    

    % delete empty station breaks (which are there if sinex is badly written)
    if deleteEmptyBreaks
        % for all available stations/epo
        for k=1:length(solEst.stat)
            if ~isempty(solEst.stat(k).break)
                allFields=fieldnames(solEst.stat(k).break); % fields of break
                nBr=length(solEst.stat(k).break);
                breaksEmpty=zeros(nBr,1);
                for kBr=1:nBr
                    nFieldsEmpty=0;
                    for kF=1:length(allFields)
                        if isempty(solEst.stat(k).break(kBr).(allFields{kF}))
                            nFieldsEmpty=nFieldsEmpty+1;
                        end
                    end
                    if nFieldsEmpty==length(allFields)
                        breaksEmpty(kBr,1)=1;
                    end
                end
                % delete the empty breaks
                solEst.stat(k).break(logical(breaksEmpty))=[];
            end
        end
    end
     
end

% --------------------------------------------------------------
function siteId=readSiteIdBlock(fid)
    curBlockName='SITE/ID';
    
    curl=fgetl(fid);
    if isempty(curl); curl=' '; end
%     kSite=0;
%     preallLength=5000;
%     code=zeros(preallLength,1); pt=cell(preallLength,1); domes=pt;
%     t=pt; station=pt; descr=pt;
    
    while ~strcmpi(deblank(curl),['-',curBlockName])
        % if not a comment line
        if ~strcmpi(curl(1),'*')
            
            % if we already have that code
            if sum(siteCodesGlobal==str2double(curl(2:5)))>0
                curCodeInd=find( siteCodesGlobal==str2double(curl(2:5)) );
            else
                curCodeInd=sum(siteCodesGlobal~=0)+1;
                siteCodesGlobal(curCodeInd)=str2double(curl(2:5)); %  save it to "global" var -> this is the "main" vector for sitecodes
                siteId(curCodeInd).code=siteCodesGlobal(curCodeInd);  
                siteId(curCodeInd).pt=curl(7:8);
                siteId(curCodeInd).domes=curl(10:18);
                siteId(curCodeInd).t=curl(20);
                %siteId(curCodeInd).station=curl(22:29);
                siteId(curCodeInd).descr=curl(22:43);
                siteId(curCodeInd).apprLon=...
                    ( str2double(curl(46:47))+str2double(curl(49:50))/60+...
                    str2double(curl(52:55))/60/60 )*pi/180;
                siteId(curCodeInd).apprLat=...
                    ( str2double(curl(58:59))+str2double(curl(61:62))/60+...
                    str2double(curl(64:67))/60/60 )*pi/180;
            end
            
%             kSite=kSite+1;
%             code(kSite)=str2double(curl(2:5));
%             pt{kSite}=curl(7:8);
%             domes{kSite}=curl(10:18);
%             t{kSite}=curl(20);
%             station{kSite}=curl(22:29);
%             descr{kSite}=curl(31:43);
            % TO DO: Approx coordinates!
            
        end
        
        % read next line
        curl=fgetl(fid);
        if isempty(curl); curl=' '; end
        
    end
    
    % delete not wanted lines
%     if preallLength>kSite
%         code(kSite+1:end)=[];
%         pt(kSite+1:end)=[];
%         domes(kSite+1:end)=[];
%         t(kSite+1:end)=[];
%         station(kSite+1:end)=[];
%         descr(kSite+1:end)=[];
%     end
% 
%     aaaa=5;
    
end

end
