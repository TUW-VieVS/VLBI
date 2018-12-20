% #########################################################################
% #     calc_separation_angle
% #########################################################################
%
% DESCRIPTION
%   Calculates the separation angle (angular distance) between two
%   target-directions defined by azimuth an elevation.
%
% CREATED  
%   2015-03-18     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%
%
% INPUT
% - az1                 - Azimuth of target 1 [rad]
% - el1                 - Elevartion of target 1 [rad]
% - az2                 - Azimuth of target 2 [rad]
% - el2                 - Elevartion of target 2 [rad]
%
%
% OUTPUT
% - sep_angle           - separation angle [rad]
%
% CHANGES:
%
%

function [sep_angle] = calc_separation_angle(az1, el1, az2, el2)

    % Angular distance between sun (v_sun) and source (v_src)
    
    v_1 = [cos(el1)*cos(az1);cos(el1)*sin(az1);sin(el1)]; % Unity vector 1
    v_2 = [cos(el2)*cos(az2);cos(el2)*sin(az2);sin(el2)]; % Unity vector 2

    sep_angle = acos((v_1'*v_2)); % / (norm(v_src) * norm(v_sun))); % Point product

return;
