function [R]=rotm(angle,type)

%#########################
% VECTORISED
%#########################

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gives the rotational matrix
%
% Input:    angle ..... [rad]
%           type  ..... type of rotation (1,2,3)
%                        1 ... x
%                        2 ... y
%                        3 ... z
%
% Output:   R     ..... rotational matrix (3,3,n)
% 
% Changelog:
% 09.05.16 M. Schartner: vectorised 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = numel(angle);
ca = cos (angle);
sa = sin (angle);

R = zeros(3,3,n);
switch type
    case 1
            R(1,1,:) = 1;
            R(2,2,:) = ca;
            R(2,3,:) = sa;
            R(3,2,:) = -sa;
            R(3,3,:) = ca;
    case 2
            R(1,1,:) = ca;
            R(1,3,:) = -sa;
            R(2,2,:) = 1;
            R(3,1,:) = sa;
            R(3,3,:) = ca;
    case 3
            R(1,1,:) = ca;
            R(1,2,:) = sa;
            R(2,1,:) = -sa;
            R(2,2,:) = ca;
            R(3,3,:) = 1;
end
