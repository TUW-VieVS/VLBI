% #########################################################################
% #     getEstimatedFields
% #########################################################################
%
% DESCRITPION
% This function returns the fields that are estimated in the sessions given
% in "indices". If only one session is plotted, this variable consists of a
% scalar (e.g. 5), otherwise it is a vector, e.g. [1,2,3,4,5,6,7,8].
%
% AUTHOR 
%   Matthias Madzak
%
% INPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%   chosenPanel  Number of the chosen panel/folder (1, 2, or 3)
%   indices      Indices of the selected sessions
%
% OUTPUT
%   estimatedFields   Estimates fields in the selected sessions
%
% CHANGES
%   2014-08-18, A. Hellerschmied: Adaptions to plot source coordinates in
%       parameter plotting window.

function estimatedFields=getEstimatedFields(handles, indices, chosenPanel)

% get fieldnames of first x_ files
fieldsOfFirstX_=fieldnames(handles.plot.x_files(chosenPanel).x_(1));
nParams=size(fieldsOfFirstX_,1);

% preallocating (one line=one session, one col=one parameter)
temp_fieldsEstimated=ones(length(indices), size(fieldsOfFirstX_,1));

% for all sessions
for k=1:length(indices)
    curInd=indices(k);

    % for all (possible estimated) parameters
    for iPara=1:nParams
        % see if the parameter was estimated
        if ~isfield(handles.plot.x_files(chosenPanel).x_(curInd).(fieldsOfFirstX_{iPara}), 'col')
            temp_fieldsEstimated(k,iPara)=0;
        elseif isempty([handles.plot.x_files(chosenPanel).x_(curInd).(fieldsOfFirstX_{iPara}).col])
            temp_fieldsEstimated(k,iPara)=0;
        end
    end
end

% get parameters that were estimated in minimum 1 sessions
estimatedFields=fieldsOfFirstX_(sum(temp_fieldsEstimated,1)>=1);