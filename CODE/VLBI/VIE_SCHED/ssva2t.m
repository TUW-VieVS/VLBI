% #########################################################################
% #     ssva2t.m
% #########################################################################
%
% DESCRITPION
%   Calculate the slew time considering acceleration and deceleration [sec]. 
%
% AUTHOR 
%   SUN Jing (2010-04-06)
%
% INPUT
%   s       angular distance [deg]
%   v       slew rate [deg/sec]
%   a       acceleration [deg/sec^2]
%
% OUTPUT
%   t       slew time [sec]
%
% COUPLING
%
%
% CHANGES
%   2015-03-02, A. Hellerschmied: Header added.


function [t] = ssva2t(s, v, a)

t1 = v/a;
s1 = 2*(a*t1*t1/2);   % acceleration + deceleration

if (s <= s1)
    t = 2*sqrt(s/a);
else
    t = 2*t1 + (s-s1)/v;
end


