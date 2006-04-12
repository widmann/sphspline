% plotsserpmapproj() - Project channel loaction coordinates to plane,
%                      helper function for pop_plotsserpmap
%
% Usage:
%   >> [F, ePrime, x, convexHull, headRadius] =
%          plotsserpmapproj(E, proj, gridScale);
%
% Inputs:
%   E             - nbchan by 3 matrix with cartesian channel location
%                   coordinates x, y, z of measured channels
%
% Optional inputs:
%   proj          - string type of projection. 'stereographic',
%                   'orthographic', 'gnomonic', 'equidistant', or
%                   'equiareal' {default 'equiareal'}
%   gridScale     - scalar integer size of interpolated scalp map data
%                   matrix {default: 61}
%
% Outputs:
%   F             - gridScale ^ 2 by 3 matrix with cartesian coordinates
%                   x, y, z of interpolated scalp map data matrix
%   ePrime        - nbchan by 2 matrix cartesian coordinates x, y of
%                   projected channel locations
%   x             - vector cartesian coordinates x (and y) of
%                   interpolated scalp map data matrix
%   convexHull    - vector indices of channels in ePrime which define
%                   the convex hull of ePrime
%   headRadius    - scalar radial distance of projected equator from
%                   vertex
%
% Author: Andreas Widmann, University of Leipzig, 2006
%
% See also:
%   pop_plotsserpmap(), plotsserpmap()

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

function [F, ePrime, x, convexHull, headRadius] = plotsserpmapproj(E, proj, gridScale)

% Arguments and defaults
if nargin < 3 || isempty(gridScale)
    gridScale = 61;
end
if nargin < 2 || isempty(proj)
    proj = 'equiareal';
end
if nargin < 1
    error('Not enough input arguments')
end
shrinkFactor = (gridScale - 1) / (gridScale - 5);

% Projection of channel locations
[th, phi] = cart2sph(E(:, 1), E(:, 2), E(:, 3)); % 3D cartesian to spherical
[ePrime(:, 1), ePrime(:, 2)] = mapproj(['sph2' proj], th, phi); % spherical to 2D cartesian
ePrime = ePrime(~sum(isnan(ePrime), 2), :); % Remove channels which could not be projected, e.g. because they are on the distant hemisphere

% Meshgrid
x = linspace(-max(abs(ePrime(:))) * shrinkFactor, max(abs(ePrime(:))) * shrinkFactor, gridScale);
[xMesh, yMesh] = meshgrid(x, x);

% Convex hull
convexHull = convhull(ePrime(:, 1), ePrime(:, 2));
inConvexHull = inpolygon(xMesh(:), yMesh(:), ePrime(convexHull, 1) * shrinkFactor, ePrime(convexHull, 2) * shrinkFactor);

% Back projection and F matrix
[th, phi] = mapproj([proj '2sph'], xMesh(inConvexHull), yMesh(inConvexHull)); % 2D cartesian to spherical
F = nan(gridScale ^ 2, 3);
[F(inConvexHull, 1) F(inConvexHull, 2) F(inConvexHull, 3)] = sph2cart(th, phi, 1); % spherical to 3D cartesian

% Head radius
headRadius = mapproj(['sph2' proj], -pi/2, 0);
