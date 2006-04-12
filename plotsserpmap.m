% plotsserpmap() - Plot map in current axis, helper function for
%                  pop_plotsserpmap
%
% Usage:
%   >> h = plotsserpmap(zData, ePrime, x, convexHull, headRadius)
%
% Inputs:
%   zData         - square matrix with interpolated data to plot
%   ePrime        - nbchan by 2 matrix cartesian coordinates x, y of
%                   projected channel locations
%   x             - vector cartesian coordinates x (and y) of of
%                   interpolated scalp map data matrix
%   convexHull    - vector indices of channels in ePrime which define
%                   the convex hull of ePrime
%   headRadius    - scalar radial distance of projected equator from
%                   vertex
%
% Outputs:
%   h             - structure with handles of the plotted objects
%
% Author: Andreas Widmann, University of Leipzig, 2006
%
% See also:
%   pop_plotsserpmap(), plotsserpmapproj()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2006 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Id$

function h = plotsserpmap(zData, ePrime, x, convexHull, headRadius)

% Arguments
if nargin < 5
    error('Not enough input arguments.');
end

% Defaults
COLOR_EEGLAB_BACKGROUND = [.93 .96 1];
ear = [[ 0.492;   0.51;  0.518; 0.5299; 0.5419;    0.54;   0.547;   0.532;    0.51;   0.484], ...
       [0.0955; 0.1175; 0.1183; 0.1146; 0.0955; -0.0055; -0.0932; -0.1313; -0.1384; -0.1199]] * 2 * headRadius; % from topoplot
nose = [[  0.09; 0.02;     0; -0.02;  -0.09], ...
        [0.4954; 0.57; 0.575;  0.57; 0.4954]] * 2 * headRadius; % from topoplot

% Figure and axis
xMin = min([ePrime(:, 1); -ear(:, 1); -headRadius]);
xMax = max([ePrime(:, 1);  ear(:, 1);  headRadius]);
yMin = min([ePrime(:, 2);             -headRadius]);
yMax = max([ePrime(:, 2); nose(:, 2);  headRadius]);
PlotBoxAspectRatio = (xMax - xMin) / (yMax - yMin);
h.axis = gca;
set(h.axis, ...
    'NextPlot', 'add', ...
    'PlotBoxAspectRatio', [PlotBoxAspectRatio 1 1], ...
    'Visible', 'off', ...
    'XLim', [xMin xMax], ...
    'YLim', [yMin yMax])

% Contour
% imagesc(x, x, zData)
h.surface = surf(x, x, zeros(size(zData)), zData, 'FaceColor', 'interp', 'LineStyle', 'none');
% [foo, h.contourf] = contourf(x, x, zData, 64, 'LineStyle', 'none');

% Contour lines
[foo, h.contour] = contour(x, x, zData, 3, 'LineColor', [0 0 0]);

% Mask convex hull
patch([x(1); x(1); x(end); x(end); x(1); ePrime(convexHull, 1)], ...
      [x(1); x(end); x(end); x(1); x(1); ePrime(convexHull, 2)], ...
      ones(length(convexHull) + 5, 1), ...
      COLOR_EEGLAB_BACKGROUND, 'EdgeColor', COLOR_EEGLAB_BACKGROUND);

% Electrodes
plot3(ePrime(:, 1), ePrime(:, 2), ones(size(ePrime, 1), 1), 'k.')

% Head, ears and nose
plot3(sin(0 : 2 * pi / 100 : 2 * pi) * headRadius, cos(0 : 2 * pi / 100 : 2 * pi) * headRadius, ones(101, 1), 'k', 'linewidth', 2)
plot3(ear(:, 1), ear(:,2), ones(size(ear, 1), 1), 'k', 'linewidth', 2)
plot3(-ear(:, 1), ear(:,2), ones(size(ear, 1), 1), 'k', 'linewidth', 2)
plot3(nose(:, 1), nose(:,2), ones(size(nose, 1), 1), 'k', 'linewidth', 2)

% drawnow
