% sserp() - Perform spherical spline interpolation
%
% Usage:
%   >> erpData = sserp(C, g, type)
%
% Inputs:
%   C         - nbchan + 1 by pnts matrix weights
%   g         - locs by nbchan matrix g(cos(E, F))
%
% Optional inputs:
%   type      - string type of interpolation. 'sp' (scalp potential,
%               µV), 'scd' (scalp current density, scaled to mA/m^3), or
%               'lap' (surface Laplacian, scaled to mV/m^2) {default:
%               'sp'}
%
% Outputs:
%   erpData   - locs by pnts matrix interpolated data
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
%   sserpgfcn(), sserpweights()

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

function erpData = sserp(C, g, type)

% Arguments
if nargin < 3 || isempty(type)
    type = 'sp';
end
if nargin < 2 || isempty(g) || isempty(C)
    error('Not enough input arguments.')
end

% Constants
CONDUCTIVITY = -0.45; % Siemens/meter
HEAD_RADIUS = 0.1; % meter
SCALE = 1e-3; % µA/m^3 to mA/m^3 and µV/m^2 to mV/m^2

% Interpolation
switch type
    case 'sp'
        % Perrin et al., 1989, eqn. (1)
        erpData = C(ones(1, size(g, 1)), :) + g * C(2:end, :);
    case 'scd'
        % Perrin et al., 1989, eqn. (5), negative 2D spherical laplacian
        erpData = -g * C(2:end, :) / HEAD_RADIUS ^ 2 * CONDUCTIVITY * SCALE;
    case 'lap'
        % Perrin et al., 1989, eqn. (5), negative 2D spherical laplacian
        erpData = -g * C(2:end, :) / HEAD_RADIUS ^ 2 * SCALE;
    otherwise
        error('Unrecognized or ambiguous interpolation type specified.');
end
