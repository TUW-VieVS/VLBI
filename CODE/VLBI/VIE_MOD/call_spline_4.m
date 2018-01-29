function y = call_spline_4(xx, yy, x)

% call_spline_4.m
% This script first cuts the available sample points to 4 values before and
% 4 values after the respective interpolant value, in order to save 
% computation time. Then, the actual spline interpolation is carried out.
% If there are no 4 points available, then just 3 are used.
%
% Coded for VieVS: 
% 27 Apr 2016 by Daniel Landskron
%
%   Revision:
%   27 Jan 2017 by Daniel Landskron: exceptions made if there are less than 4 values before/after, also extrapolation
%
% *************************************************************************



num_points = length(x);
y = zeros(num_points,1);

for i_point = 1:num_points
    
    % new selection of proper intervals (4 before , 4 after -> 8 values always)
    indSmaller = flipud(find(xx<=x(i_point)));
    indLarger = find(xx>=x(i_point));
    if length(indSmaller)>=4
        indFirst = indSmaller(4);
    elseif length(indSmaller)==3
        indFirst = indSmaller(3);
    elseif length(indSmaller)==2
        indFirst = indSmaller(2);
    elseif length(indSmaller)==1
        indFirst = indSmaller(1);
    else
        disp('too few sample points before the interpolant to perform a spline interpolation, extrapolated instead');
        indFirst = indLarger(1);
    end
    if length(indLarger)>=4
        indLast = indLarger(4);
    elseif length(indLarger)==3
        indLast = indLarger(3);
    elseif length(indLarger)==2
        indLast = indLarger(2);
    elseif length(indLarger)==1
        indLast = indLarger(1);
    else
        disp('too few sample points after the interpolant to perform a spline interpolation, extrapolated instead');
        indLast = indSmaller(1);
    end
    
    xx_temp = xx(indFirst:indLast);
    yy_temp = yy(indFirst:indLast);
    
    % do the spline interpolation for every value x
    y(i_point) = spline(xx_temp, yy_temp, x(i_point));
    
end