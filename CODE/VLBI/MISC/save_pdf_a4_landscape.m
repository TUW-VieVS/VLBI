% #########################################################################
% #     save_pdf_a4_landscape
% #########################################################################
% DESCRIPTION
%   Function to save plots as PDF file. Different MATLAB releases are considered.
% Paper format: A4, landscape, fill page
%
% CREATED  
%    A. Hellerschmied, 2016-11-09
%
% COUPLING
%
% INPUT
% - h_fig               - figure handle
% - source              - path of output file (ralative to work dor., or absolut), "/" in the end - string
% - filename            - Name of PDF file - string
%
%
% OUTPUT
%
function save_pdf_a4_landscape(h_fig, filepath, filename)

if strcmp(filepath(1), '.')
    tmp_str = '';
else
    tmp_str = '.';
end

% Check, if the outpur dir. exists:
if ~exist( [tmp_str, filepath], 'dir')
   mkdir([tmp_str, filepath]) 
end

mallab_release_str = version('-release');
release = str2double(mallab_release_str(1:4));
if release >= 2016
    % From MATLAB version 2016a
    h_fig.PaperOrientation = 'landscape';
    print(h_fig,'-fillpage', '-dpdf','-r600',  [tmp_str, filepath, filename]);
else
    % Older Matlab versions:
    h_fig.PaperOrientation = 'landscape';
    h_fig.PaperUnits = 'normalized';
    h_fig.PaperPosition = [0 0 1 1];

    print(h_fig,'-dpdf','-r600', [tmp_str, filepath, filename]);
end
    
return