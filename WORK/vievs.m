% #########################################################################
% #     vievs
% #########################################################################
%
% DESCRIPTION
%   This function is used to start the VieVS software.
%
%   Uusage/input parameter options:
%   - vievs                 : Runs the latest version of vievs (both setup and data processing part)
%   - vievs('SETUP')        : Runs the setup part of vievs
%   - vievs('BATCH')        : Runs the data processing part of vievs
%
% CREATED
%   ????
%
% REFERENCES
%
%
% COUPLING
%   - vievs_setup
%   - vie_batch
%
%
% INPUT
% => Find information in the DESCRIPTION section above!
%
% OUTPUT
%
% CHANGES
% - 2015-12-18, A. Hellerschmied:   - Updated for VieVS 3.0
%                                   - Downward compatibility to VieVS 1.x not any longer provided
%                                   - Error in the initial removing of existing .../COMPILE/.. directories fixed.
% - 2016-10-10, A. Girdiuk:		    - Linux compatibility added for path
% - 2017-05-31, A. Hellerschmied:   - Bug fixed at removing COMPILE directories from the MATLAB search patch
% - 2017-05-31, A. Hellerschmied:   - Updated for VieVS 3.1
% - 2018-01-18, A. Hellerschmied    - Revised for handling of VieVS with GIT
%
function vievs(varargin)

% ##### Init.: #####
close all;
flag_runsetup = 0;
flag_runbatch = 0;
clc;

% ##### Add path of the /COMMON/ and /CODE/ directories to the Matlab
% search path ####
addpath(genpath('../../COMMON/'));
addpath(genpath('../CODE/'));

% ##### Give a warning if the COMMON directory does not exist #####
if ~exist('../../COMMON/', 'dir')
    msgbox('ERROR: The /COMMON directory does not exist - please clone/copy the COMMON repository to the level of your VLBI repository!!!');
end

% ##### Check if the VLBI_OPT directory exists, if not create it ####
OptFolder='../../VLBI_OPT/';
if ~exist(OptFolder, 'dir')
    mkdir(OptFolder);
    mkdir([OptFolder,'DEFAULT/']);
end

% ##### Check optional input arguments #####
for a = 1 : nargin
    if ischar(varargin{a})
        if strcmpi(varargin{a},'SETUP')
            flag_runsetup = 1;
        elseif  strcmpi(varargin{a},'BATCH')
            flag_runbatch = 1;
        end
    end
end

if (flag_runsetup == 0) && (flag_runbatch == 0)
    flag_runsetup = 1;
    flag_runbatch = 0;
end

% ##### GUI mode (run VIE_SETUP) #####
if flag_runsetup
    vie_setup
end

% ##### Batch mode (run vie_batch only) #####
if flag_runbatch
    vie_batch
end

