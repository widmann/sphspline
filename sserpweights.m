% sserpweights() - Compute C matrix used for spherical spline
%                  interpolation 
%
% Usage:
%   >> C = sserpweights(data, Ginv)
%
% Inputs:
%   data  - nbchan by pnts matrix of measured data
%   Ginv  - matrix padded inverse of g(cos(E, E))
%
% Outputs:
%   C     - nbchan + 1 by pnts matrix weights
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
%   sserpgfcn(), sserp()

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

function C = sserpweights(data, Ginv)

% Arguments
if nargin < 2
    error('Not enough input arguments.')
end

% C matrix
C = Ginv * [zeros(1, size(data, 2)); double(data)]; % Perrin et al., 1989, eqn. (2)
