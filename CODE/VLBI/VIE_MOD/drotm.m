%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gives the derivative of the rotational matrix
%
% Input:    angle ..... [rad]
%           type  ..... type of rotation (1,2,3)
%                        1 ... x
%                        2 ... y
%                        3 ... z
%
% Output:   R     ..... rotational matrix (3,3,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%#########################
% VECTORISED
%#########################
%   Revision:
%	08 Aug 2016 by A. Girdiuk: a loop removed
%
function [R]=drotm(angle,type)

n = length(angle);
R = zeros(3,3,n);

ca = cos (angle);
sa = sin (angle);

switch type
    case 1
            R(1,1,:) = 0;
            R(2,2,:) = -sa;
            R(2,3,:) =  ca;
            R(3,2,:) = -ca;
            R(3,3,:) = -sa;
    case 2
            R(1,1,:) = -sa;
            R(1,3,:) = -ca;
            R(2,2,:) =   0;
            R(3,1,:) =  ca;
            R(3,3,:) = -sa;
    case 3
            R(1,1,:) = -sa;
            R(1,2,:) =  ca;
            R(2,1,:) = -ca;
            R(2,2,:) = -sa;
            R(3,3,:) =   0;
end
