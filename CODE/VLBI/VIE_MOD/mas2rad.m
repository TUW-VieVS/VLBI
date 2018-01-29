function rad=mas2rad(mas)
% Converts seconds of arc to radians. Vectorized.
% Useage:  rad=as2rad(as)
% Input:   as - vector of angles in seconds of arc
% Output:  rad - vector of angles in radians

rad=(mas/1000).*pi./180./3600;
