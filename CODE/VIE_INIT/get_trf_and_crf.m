% #########################################################################
% #     get_trf_and_crf
% ########################################################################
%
% DESCRITPION
% This function reads files TRF data (suoperstation file and/or manual TRF file)
% an CRF data (supersource file).
%
%
%
% AUTHOR
%   A. Hellerschmied, 2017-02-22
%
% INPUT
%  - parameter              - VieVS parameter structure (GUI parameters, etc.)
%
% OUTPUT
%  - trf                    - TRF data structure
%  - crf                    - CRF data structure
%  - parameter              - VieVS parameter structure (GUI parameters, etc.)
%
% COUPLING
%  -
%
% CHANGES
%  - YYYY-MM-DD, <name>: <description>
%
function [trf, crf, parameter] = get_trf_and_crf(parameter)

trffile = parameter.vie_init.trf;
crffile = parameter.vie_init.crf;


% ##### TRF #####

% (1) TRF: load superstation file
% if a superstation file is given (==~manual textfile is given) - load that
if strcmp(trffile{1}(end-3:end), '.mat')
    trf = load(trffile{1});
    nam = fieldnames(trf);
    trf = eval(['trf.' nam{1}]);
else % a manual trf file is given -> BUT load default superstat file in any case!!
    trf = load('../TRF/superstation.mat');
    nam = fieldnames(trf);
    trf = eval(['trf.' nam{1}]);
end

% (2) TRF: load manual TRF (if chosen)
if strcmpi(trffile{1}(end-3:end), '.txt')
    % put in superstation file (variable trf)
    trf = manualTrfToSuperstatTrf(trf,trffile{1});
    
    % write "trfname" to trffil{2}
    trffile{2} = 'manualTrf';
    parameter.vie_init.trf{2} = trffile{2};
end


% ##### CRF #####

% (1) CRF: load supersource file
if strcmp(crffile{1}(end-3:end), '.mat')
    crf = load(crffile{1});
    nam = fieldnames(crf);
    crf = eval(['crf.', nam{1}]);
else
    error('Manual CRF not supported yet!')
end