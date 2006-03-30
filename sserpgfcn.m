% sserpgfcn() - Compute G matrices used for spherical spline
%               interpolation
%
% Usage:
%   >> [Ginv, g, G] = sserpgfcn(E, F, type, lambda, nTerms, m)
%
% Inputs:
%   E         - nbchan by 3 matrix with cartesian channel location
%               coordinates x, y, z of measured channels
%   F         - locs by 3 matrix with cartesian channel location
%               coordinates x, y, z to interpolate in columns
%
% Optional inputs:
%   type      - string type of interpolation. 'sp' (scalp potential),
%               'scd' (scalp current density), or 'lap' (surface
%               Laplacian) {default 'sp'}
%   lambda    - scalar smoothing factor (commonly used values are 1e-7
%               for sp and 1e-5 for scd) {default 0}
%   nTerms    - scalar int > 0 number of terms {default 50}
%   m         - scalar int > 1 m {default 4}
%
% Outputs:
%   Ginv      - nbchan + 1 by nbchan + 1 matrix padded inverse of
%               g(cos(E, E)) 
%   g         - locs by nbchan matrix g(cos(E, F))
%   G         - nbchan by nbchan matrix g(cos(E, E))
%
% Note:
%   Recursive version of gfcn is tens of times faster than using MATLAB
%   legendre function. Sserpgfcn implicitly understands that channel
%   coordinates in E and F are located on the surface of a unit sphere.
%   It will return wrong results if they are not. Use unitsph() to
%   project channel coordinates to a unit sphere surface.
%
% References:
%   [1] Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F.
%       (1989). Spherical splines for scalp potential and current
%       density mapping. Electroencephalography and Clinical
%       Neurophysiology, 72, 184-187
%   [2] Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F.
%       (1990). Corrigenda EEG 02274. Electroencephalography and
%       Clinical Neurophysiology, 76, 565
%   [3] Kayser, J., & Tenke, C. E. (2006). Principal components analysis
%       of Laplacian waveforms as a generic method for identifying ERP
%       generator patterns: I. Evaluation with auditory oddball tasks.
%       Clinical Neurophysiology, 117, 348-368
%   [4] Weber, D. L. (2001). Scalp current density and source current
%       modelling. Retrieved March 26, 2006, from
%       dnl.ucsf.edu/users/dweber/dweber_docs/eeg_scd.html
%   [5] Ferree, T. C. (2000). Spline Interpolation of the Scalp EEG.
%       Retrieved March 26, 2006, from
%       www.egi.com/Technotes/SplineInterpolation.pdf 
%   [6] Ferree, T. C., & Srinivasan, R. (2000). Theory and Calculation
%       of the Scalp Surface Laplacian. Retrieved March 26, 2006, from
%       http://www.egi.com/Technotes/SurfaceLaplacian.pdf
%
% Author: Andreas Widmann, University of Leipzig, 2006
%
% See also:
%   sserpweights(), sserp(), unitsph()

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

function [Ginv, g, G] = sserpgfcn(E, F, type, lambda, nTerms, m)

if nargin < 6 || isempty(m)
    m = 4;
end
if nargin < 5 || isempty(nTerms)
    nTerms = 50;
end
if nargin < 4 || isempty(lambda)
    lambda = 0;
end
if nargin < 3 || isempty(type)
    type = 'sp';
end
if nargin < 2
    error('Not enough input arguments.')
end

% Cosines, quaternion based analog to Perrin et al., 1989, eqn. (4)
x = E * E';

% G matrix
G = gfcn(x, nTerms, m);

% Pad, add smoothing constant to diagonale, and invert G
Ginv = inv([0 ones(1, size(G, 2)); ones(size(G, 1), 1) G + eye(size(G)) * lambda]);

% Cosines, quaternion based analog to Perrin et al., 1989, eqn. (4)
x = F * E';

% g matrix
switch type
    case 'sp'
        g = gfcn(x, nTerms, m); % Perrin et al., 1989, eqn. (3)
    case {'scd' 'lap'}
        g = gfcn(x, nTerms, m - 1); % Perrin et al., 1990, eqn. after eqn. (5)
    otherwise
        error('Unrecognized or ambiguous interpolation type specified.');
end

function [G] = gfcn(x, nTerms, m)
    P = cat(3, ones(size(x)), x);
    G = 3 / 2 ^ m * P(:, :, 2);
    for n = 2:nTerms
        % nth degree legendre polynomial; Perrin et al., 1989, eqn. after eqn. (4)
        P(:, :, 3) = ((2 * n - 1) * x .* P(:, :, 2) - (n - 1) * P(:, :, 1)) / n;
        P(:, :, 1) = [];
        % Perrin et al., 1989, eqn. (3)
        G = G + (2 * n + 1) / (n ^ m * (n + 1) ^ m) * P(:, :, 2);
    end
    G = G / (4 * pi);
