% -------------------------------------------------------------------------
%
%                              skyplot.m
%
%   Creates a skyplot for locations gives in AZ/EL.
%
%   As basis the following MATLAB function was used:
%   http://read.pudn.com/downloads120/sourcecode/others/511876/sofGPS/include/skyPlot.m__.htm
%
%   
%
%
%   Author: 
%       Andreas Hellerschmied (heller182@gmx.at), 25.02.2013
%   
%   changes       :
%   - 2014-08-28: A. Hellerschmied: return axis-handle (hAxis); plot
%       markers as "timescale" in skyplot. 
%   - 2015-11-26: A. Hellerschmied: Return handles "h_marker, "h_sat_num_txt" and "h_sat_marker" as argument
%   - 2016-04-14: A. Hellerschmied: Bug-fix: line_style parameter was missing for nargin = 4 
%           
%
%   inputs        :
%   - az            :   contains satellite azimuth angles. It is a 2D 
%                       matrix. One line contains data of one satellite.
%                       The columns are the calculated azimuth values. 
%   - el            :   contains satellite elevation angles. It is a 2D 
%                       matrix. One line contains data of one satellite. 
%                       The columns are the calculated elevations.
%   - prn           :   a row vector containing PRN numbers of the 
%                       satellites.
%   - marker_int    :   Interval the "time-scale" markers (measured in 
%                       initial propagation intervals).
%   - line_style    :   line style of the plot. The same style will be 
%                       used to plot all satellite positions (including 
%                       color).
%     
%
%   outputs       :
%   - hpol          : plot handles
%   - hAxis         : axis handles
%   - plotted_sats  : logical vector (1 and 0), indictaiong which
%                       satellites have been plotted (are visible).
%   - h_marker      : timescale marker handles (array)
%   - h_sat_num_txt         : satellite number text handles (array)
%   - h_sat_marker          : satellite marker handles (array)
%    
%
%   locals        :
% 
%
%   coupling      :
%   
%
%   references    : Modified from function referenced below!
%
%-------------------------------------------------------------------------


function [hpol, hAxis, plotted_sats, h_marker, h_sat_num_txt, h_sat_marker] = skyplot(varargin) 

%Function plots "sky view" from the receiver perspective.  
% 
%h = skyPlot(AZ, EL, PRN, line_style) 
% 
%   Inputs: 
%       AZ              - contains satellite azimuth angles. It is a 2D 
%                       matrix. One line contains data of one satellite.
%                       The columns are the calculated azimuth values. 
%       EL              - contains satellite elevation angles. It is a 2D 
%                       matrix. One line contains data of one satellite. 
%                       The columns are the calculated elevations. 
%       PRN             - a row vector containing PRN numbers of the 
%                       satellites. 
%       line_style      - line style of the plot. The same style will be 
%                       used to plot all satellite positions (including 
%                       color).  
%   Outputs: 
%       h               - handle to the plot 
%-------------------------------------------------------------------------- 
%                           SoftGNSS v3.0 
%  
% Copyright (C) Darius Plausinaitis and Kristin Larson 
% Written by Darius Plausinaitis and Kristin Larson 
%-------------------------------------------------------------------------- 
%This program is free software; you can redistribute it and/or 
%modify it under the terms of the GNU General Public License 
%as published by the Free Software Foundation; either version 2 
%of the License, or (at your option) any later version. 
% 
%This program is distributed in the hope that it will be useful, 
%but WITHOUT ANY WARRANTY; without even the implied warranty of 
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
%GNU General Public License for more details. 
% 
%You should have received a copy of the GNU General Public License 
%along with this program; if not, write to the Free Software 
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%
%USA. 
%-------------------------------------------------------------------------- 
 
%CVS record: 
%$Id: skyPlot.m,v 1.1.2.5 2006/08/18 11:41:57 dpl Exp $ 
 
%% Check arguments and sort them ========================================== 
[hAxis, args, nargs] = axescheck(varargin{:}); 
 
if nargs < 3 || nargs > 5
    error('Requires 3 to 5 data arguments.') 
elseif nargs == 3 
    [az, el, prn]   = deal(args{1:3}); 
    line_style      = 'auto';   
    marker_int      = 0;
elseif nargs == 4 
    [az, el, prn, marker_int] = deal(args{1:4});
    line_style      = 'auto';   
else % nargs == 5
    [az, el, prn, marker_int, line_style] = deal(args{1:5});
end 
 
if ischar(az) || ischar(el) || ischar(prn) 
    error('AZ and EL must be numeric.'); 
end 
 
if ~isequal(size(az), size(el)) 
    error('AZ and EL must be same size.'); 
end 
 
%% Prepare axis ===========================================================
hAxis = newplot(hAxis); 
 
%--- Get x-axis text color so grid is in same color ----------------------- 
tc = get(hAxis, 'xcolor'); 
 
hold(hAxis, 'on'); 
 
%--- Plot white background ------------------------------------------------ 
rectangle('position', [-90, -90, 180, 180], ... 
          'Curvature', [1 1], ... 
          'facecolor', 'white', ... 
          'edgecolor', tc); 
 
%% Plot spokes ============================================================ 
 
%--- Find spoke angles ---------------------------------------------------- 
% Only 6 lines are needed to divide circle into 12 parts 
th = (1:6) * 2*pi / 12; 
 
%--- Convert spoke end point coordinate to Cartesian system --------------- 
cst = cos(th); snt = sin(th); 
cs = [cst; -cst]; 
sn = [snt; -snt]; 
 
%--- Plot the spoke lines ------------------------------------------------- 
line(90*sn, 90*cs, 'linestyle', ':', 'color', tc, 'linewidth', 0.5, ... 
    'handlevisibility', 'off'); 
 
%% Annotate spokes in degrees ============================================= 
rt = 1.1 * 90; 
 
for i = 1:max(size(th)) 
 
    %--- Write text in the first half of the plot ------------------------- 
    text(rt*snt(i), rt*cst(i), int2str(i*30), ... 
        'horizontalalignment', 'center', 'handlevisibility', 'off'); 
 
    if i == max(size(th)) 
        loc = int2str(0); 
    else 
        loc = int2str(180 + i*30); 
    end 
 
    %--- Write text in the opposite half of the plot ---------------------- 
    text(-rt*snt(i), -rt*cst(i), loc, ... 
        'handlevisibility', 'off', 'horizontalalignment', 'center'); 
end 
 
%% Plot elevation grid ==================================================== 
 
%--- Define a "unit" radius circle ---------------------------------------- 
th = 0 : pi/50 : 2*pi; 
xunit = cos(th); 
yunit = sin(th); 
 
%--- Plot elevation grid lines and tick text ------------------------------ 
for elevation = 0 : 15 : 90 
    elevationSpherical = 90*cos((pi/180) * elevation); 
 
    line(yunit * elevationSpherical, xunit * elevationSpherical, ... 
        'lineStyle', ':', 'color', tc, 'linewidth', 0.5, ... 
        'handlevisibility', 'off'); 
 
    text(0, elevationSpherical, num2str(elevation), ... 
        'BackgroundColor', 'white', 'horizontalalignment','center', ... 
        'handlevisibility', 'off'); 
end 
 
%--- Set view to 2-D ------------------------------------------------------ 
view(0, 90); 
 
%--- Set axis limits ------------------------------------------------------ 
%save some space for the title 
axis([-95 95 -90 101]); 
 
%% Transform elevation angle to a distance to the center of the plot ------ 
elSpherical = 90*cos(el * pi/180); 
 
%--- Transform data to Cartesian coordinates ------------------------------ 
yy = elSpherical .* cos(az * pi/180); 
xx = elSpherical .* sin(az * pi/180); 
 
%% Plot data on top of the grid =========================================== 
 
if strcmp(line_style, 'auto') 
    %--- Plot with "default" line style ----------------------------------- 
    hpol = plot(hAxis, xx', yy', '.-'); 
else 
    %--- Plot with user specified line style ------------------------------ 
    % The same line style and color will be used for all satellites 
    hpol = plot(hAxis, xx', yy', line_style, 'linewidth', 2); 
end 

%--- Transpose xx and yy vectors for easier further use -------------------
xx = xx';
yy = yy';
 
%--- Mark the last position of the satellite and Place PRN numbers --------
plotted_sats = zeros(length(yy(1,:)),1);
for i_sat = 1 : length(yy(1,:))
    % Select datapoints for markers:
    xx_temp = xx(~isnan(xx(:,i_sat)), i_sat);
    yy_temp = yy(~isnan(yy(:,i_sat)), i_sat);
    % Plot markers in the same color as the related lines:
    if ~isempty(xx_temp) || ~isempty(yy_temp)
        plotted_sats(i_sat) = 1;
        h_sat_marker(i_sat) = plot(hAxis, xx_temp(end,:), yy_temp(end,:), 'o', 'Color', get(hpol(i_sat), 'Color'), 'MarkerSize', 7);

        % Place satellite PRN numbers at the latest position, if available:
        if (length(prn) ==  length(yy(1,:)))
             h_sat_num_txt(i_sat) = text(xx_temp(end,:), yy_temp(end,:), ['  ', int2str(prn(i_sat))], 'Color', get(hpol(i_sat), 'Color'));
        end
    end
end

%--- Plot "timescale markers"  --------------------------------------------
if marker_int > 0
    for i_sat = 1 : length(yy(1,:))
        if plotted_sats(i_sat) == 1
            % Select datapoints for markers:
            xx_marker = xx((~isnan(xx(:,i_sat)) & (mod(find(xx(:,i_sat)), marker_int) == 0)), i_sat);
            yy_marker = yy((~isnan(yy(:,i_sat)) & (mod(find(yy(:,i_sat)), marker_int) == 0)), i_sat);
            % Plot markers in the same color as the related lines:
            if ~isempty(xx_marker) || ~isempty(yy_marker)
                h_marker(i_sat) = plot(hAxis, xx_marker, yy_marker, '*', 'Color', get(hpol(i_sat), 'Color'));
            end
        end
    end
end

if ~exist('h_marker', 'var') 
    h_marker = []; 
end
if ~exist('h_sat_num_txt', 'var') 
    h_sat_num_txt = []; 
end
if ~exist('h_sat_marker', 'var') 
    h_sat_marker = []; 
end

%--- Make sure both axis have the same data aspect ratio ------------------ 
axis(hAxis, 'equal'); 
 
%--- Switch off the standard Cartesian axis ------------------------------- 
axis(hAxis, 'off'); 
