% mapproj() - Spherical projection and back projection functions
%
% Usage:
%   >> [x, y] = mapproj(type, th, phi); % projection
%   >> [th, phi] = mapproj(type, x, y); % back projection
%   >> [x, y] = mapproj(type, th, phi, th0, phi1);
%   >> [th, phi] = mapproj(type, x, y, th0, phi1);
%
% Inputs:
%   type  - char array type of projection. 'sph2stereographic',
%           'sph2orthographic', 'sph2gnomonic', 'sph2equidistant', or
%           'sph2equiareal' for projection. 'stereographic2sph',
%           'orthographic2sph', 'gnomonic2sph', 'equidistant2sph', or
%           'equiareal2sph' for  back projection.
%   th    - scalar or vector longitude in radians
%   phi   - scalar or vector latitude in radians
%   x     - scalar or vector cartesian x
%   y     - scalar or vector cartesian y
%
% Optional inputs:
%   th0   - scalar central longitude {default pi}
%   phi1  - scalar central latitude {default pi / 2}
%
% Outputs:
%   x     - scalar or vector projected cartesian x
%   y     - scalar or vector projected cartesian y
%   th    - scalar or vector back projected longitude in radians
%   phi   - scalar or vector back projected latitude in radians
%
% References:
%   [1] Weisstein, E. W. (1999). Stereographic Projection. Retrieved
%       March 19, 2006, from
%       http://mathworld.wolfram.com/StereographicProjection.html
%   [2] Weisstein, E. W. (1999). Orthographic Projection. Retrieved
%       March 19, 2006, from
%       http://mathworld.wolfram.com/OrthographicProjection.html
%   [3] Weisstein, E. W. (1999). Gnomonic Projection. Retrieved March
%       19, 2006, from
%       http://mathworld.wolfram.com/GnomonicProjection.html
%   [4] Weisstein, E. W. (1999). Azimuthal Equidistant Projection.
%       Retrieved March 19, 2006, from
%       http://mathworld.wolfram.com/AzimuthalEquidistantProjection.html
%   [5] Weisstein, E. W. (1999). Lambert Azimuthal Equal-Area
%       Projection. Retrieved  March 19, 2006, from
%       http://mathworld.wolfram.com/LambertAzimuthalEqual-
%       AreaProjection.html
%
% Author: Andreas Widmann, University of Leipzig, 2006

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

function varargout = mapproj(proj, varargin)

    % Arguments and defaults
    if nargin < 3
        error('Not enough input arguments.')
    end
    if nargin < 4
        varargin{end + 1} = pi; % Nose on top
    end
    if nargin < 5
        varargin{end + 1} = pi / 2; % Vertex in center
    end

    proj = str2func(proj);
    [varargout{1:nargout}] = proj(varargin{:});

% Projection
function [x, y] = sph2stereographic(th, phi, th0, phi1)
    k = 2 ./ (1 + sin(phi1) * sin(phi) + cos(phi1) * cos(phi) .* cos(th - th0));
    [x, y] = sph2proj(k, th, phi, th0, phi1);

function [x, y] = sph2orthographic(th, phi, th0, phi1)
    phi(phi < 0) = NaN;
    k = 1;
    [x, y] = sph2proj(k, th, phi, th0, phi1);

function [x, y] = sph2gnomonic(th, phi, th0, phi1)
    phi(phi <= 0) = NaN;
    k = 1 ./ (sin(phi1) * sin(phi) + cos(phi1) * cos(phi) .* cos(th - th0));
    [x, y] = sph2proj(k, th, phi, th0, phi1);

function [x, y] = sph2equidistant(th, phi, th0, phi1)
    c = acos(sin(phi1) * sin(phi) + cos(phi1) * cos(phi) .* cos(th - th0));
    temp = c ~= 0;
    k = zeros(size(c));
    k(temp) = c(temp) ./ sin(c(temp));
    [x, y] = sph2proj(k, th, phi, th0, phi1);

function [x, y] = sph2equiareal(th, phi, th0, phi1)
    k = sqrt(2 ./ (1 + sin(phi1) * sin(phi) + cos(phi1) * cos(phi) .* cos(th - th0)));
    [x, y] = sph2proj(k, th, phi, th0, phi1);

% Projection main function
function [x, y] = sph2proj(k, th, phi, th0, phi1)
    x = k .* cos(phi) .* sin(th - th0);
    y = k .* (cos(phi1) * sin(phi) - sin(phi1) * cos(phi) .* cos(th - th0));

% Back projection
function [th, phi] = stereographic2sph(x, y, th0, phi1)
    rho = sqrt(x .^ 2 + y .^ 2);
    c = 2 * atan2(rho, 2);
    [th, phi] = proj2sph(rho, c, x, y, th0, phi1);

function [th, phi] = orthographic2sph(x, y, th0, phi1)
    rho = sqrt(x .^ 2 + y .^ 2);
    rho(rho > 1) = NaN;
    c = asin(rho);
    [th, phi] = proj2sph(rho, c, x, y, th0, phi1);

function [th, phi] = gnomonic2sph(x, y, th0, phi1)
    rho = sqrt(x .^ 2 + y .^ 2);
    c = atan(rho);
    [th, phi] = proj2sph(rho, c, x, y, th0, phi1);

function [th, phi] = equidistant2sph(x, y, th0, phi1)
    rho = sqrt(x .^ 2 + y .^ 2);
    c = rho;
    [th, phi] = proj2sph(rho, c, x, y, th0, phi1);

function [th, phi] = equiareal2sph(x, y, th0, phi1)
    rho = sqrt(x .^ 2 + y .^ 2);
    c = 2 * asin(rho / 2);
    [th, phi] = proj2sph(rho, c, x, y, th0, phi1);

% Back projection main function
function [th, phi] = proj2sph(rho, c, x, y, th0, phi1)
    temp = rho ~= 0;
    phi = phi1(ones(size(rho)));
    phi(temp) = asin(cos(c(temp)) .* sin(phi1) + y(temp) .* sin(c(temp)) .* cos(phi1) ./ rho(temp));
    th = th0 + atan2((x .* sin(c)), (rho * cos(phi1) .* cos(c) - y * sin(phi1) .* sin(c)));
