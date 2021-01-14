function [fb,in,ambcorr,ioncorr,wrapper_v,wrapper_k] = vgosdbinGUI2vieinit(a,b,c,d,e)
% convert gui values to vie_init values
% input:
%   a ... parameter.vie_init.vgosDb_observation_parameter
%   b ... parameter.vie_init.vgosDb_institute
%   c ... parameter.vie_init.ambiguity_correction
%   d ... parameter.vie_init.iono_correction
%   e ... parameter.vie_init.vgosDb_wrapper_version

fb=[];
in=[];
ambcorr=[];
ioncorr=[];
wrapper_v=[];
wrapper_k=[]; %% not supported yet

% vgosDb_observation_parameter
if ~isempty(a)
    if ischar(a)
        fb = {a};
    end
end

% vgosDb_institute
if ~isempty(b)
    if ischar(b)
        in = strsplit(b);
    end
end

% ambiguity_correction
if ~isempty(c)
    if c == 0
        ambcorr = 'off';
    elseif c==1
        ambcorr = 'on';
    end
end

% iono_correction
if ~isempty(d)
    if d == 0
        ioncorr = 'off';
    elseif d==1
        ioncorr = 'on';
    end
end

% vgosDb_wrapper_version
if ~isempty(e)
    wrapper_v=str2num(e);
end

end

