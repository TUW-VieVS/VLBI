% Calculate correction of post-seismic deformation
%
% CREATED: Matthias Madzak, 2016-05-10
%
% INPUT
%  mjd     modified julian date as vector (any dimension)
%  antenna vievs antenna file: must contain 'psd' field. If empty: Not psd
%          correction applied (return values are zero)
%          can be superstation file -> then, varargin must be set
%  varargin
%  antName cellst(or char array) defining the stations in superstation file
%          if that is set, 'antenna' must be superstat file (or path to mat
%          file)
%  trfName eg 'itrf2014' (only if superstation file is given in parameter
%          'antenna'
%
% OUTPUT
%  cpsd    matrix containing post-seismic deformations for all stations in
%          antenna file for all mjds
%          [3 x nScans x nStat] ... x/y/z
% CHANGES
%   2017-05-30, A. Girdiuk: check for empty array is added because mjd is before brake
%                           
function cpsd=cPostSeismDeform(mjd,antenna,varargin)

mjd=mjd(:); % make column vector

sstatFileGiven=0;
if nargin>2
    if ~isempty(varargin{1})
        sstatFileGiven=1;
        antName=cellstr(varargin{1});
        trfName=varargin{2};
        nStat=length(antName);
        if ischar(antenna)
            temps=load(antenna); f=fieldnames(temps);
            antenna=temps.(f{1}); % This is now the superstation file
        end
    end
end


if sstatFileGiven==0
    nStat=length(antenna);
end
nEp=length(mjd);

cpsd=zeros(3,nEp,nStat);

for kStat=1:nStat
    if sstatFileGiven
        % find antenna index
        foundStatLog=strcmpi({antenna.name},antName{kStat}) | ...
            strcmpi({antenna.domes},antName{kStat}) | ...
            strcmpi({antenna.CDP},antName{kStat});
        if sum(foundStatLog)==0
            fprintf('ERROR: Station not found!!\nfunction cPostSeismDeform.m\n');
        end
        if ~isfield(antenna(foundStatLog).(trfName),'psd')
            psdField=[];
        else
            psdField=antenna(foundStatLog).(trfName).psd;
        end
    else % vievs-antenna file is given as antenna
        if isfield(antenna,'psd')
            psdField=antenna(kStat).psd;
        else
            psdField = [];
        end
    end
    if ~isempty(psdField)
        if sstatFileGiven
            curx=antenna(foundStatLog).(trfName).break(1).x; %  that's a bit wrong -> could be different break!!
            cury=antenna(foundStatLog).(trfName).break(1).y;
            curz=antenna(foundStatLog).(trfName).break(1).z;
        else
            curx=antenna(kStat).x;
            cury=antenna(kStat).y;
            curz=antenna(kStat).z;
        end
        % for all breaks of psd
        for kBreak=1:length(psdField)
            mjdq=psdField(kBreak).epoch;
            curEpLog=mjd>mjdq; % epoch-logicals of current psd break
            dtq=(mjd(curEpLog)-mjdq)./365.25; % time since break [years]
            psdCorrEOfCurBreak=psdCalc(dtq,psdField(kBreak).e);
            psdCorrNOfCurBreak=psdCalc(dtq,psdField(kBreak).n);
            psdCorrUOfCurBreak=psdCalc(dtq,psdField(kBreak).u);
            
            % dENU -> dXYZ
            if ~isempty([psdCorrUOfCurBreak,psdCorrEOfCurBreak,psdCorrNOfCurBreak])
                [mjdxyz]=call_ren2xyz([mjd(curEpLog),psdCorrUOfCurBreak,...
                    psdCorrEOfCurBreak,psdCorrNOfCurBreak],...
                    [curx,cury,curz]);
            % either add to previous corr (or new corr after break)
            % I THINK -> ADD!!
                cpsd(:,curEpLog,kStat)=cpsd(:,curEpLog,kStat)+...
                mjdxyz(:,2:4)';
            end        
        end
    end
    
end

   % This function calculates the actual correction values based on dtq
    % dtq = (mjd - mjdq) / 365.25d0; where mjdq is the break epoch
    % and corrVals (ie lines of ITRF psd)
    % corrVals   eg: [2   29.69  2.4515] or [0]
    % OUTPUT: correction (usually meters)
    function corr=psdCalc(dtq,corrVals)

        corr=zeros(size(dtq));
        switch corrVals(1)
            case 0
                corr=zeros(size(dtq));
            case 1
                corr=corrVals(2).*log(1+dtq./corrVals(3));
            case 2
                corr=corrVals(2).*(1-exp(-dtq./corrVals(3)));
            case 3
                corr=corrVals(2).*log(1+dtq./corrVals(3)) + ...
                    corrVals(4).*(1-exp(-dtq./corrVals(5)));
            case 4
                corr=corrVals(2).*(1-exp(-dtq./corrVals(3))) + ...
                    corrVals(4).*(1-exp(-dtq./corrVals(5)));
        end

        corr=corr./1000; 
    
    end

end
