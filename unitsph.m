% unitsph() - Re-center channel location coordinates and project to
%             unit sphere surface
%
% Usage:
%   >> E = unitsph(E, thresh)
%
% Inputs:
%   E     - nbchan by 3 matrix with cartesian channel location
%           coordinates x, y, z
%
% Optional inputs:
%   thres - scalar threshold < abs(radius - 1) {default 1e-6}
%
% Outputs:
%   E     - nbchan by 3 matrix with cartesian channel location
%           coordinates x, y, z re-centered to best fitting sphere and
%           projected to unit sphere surface 
%
% Author: Andreas Widmann, University of Leipzig, 2006
%
% See also:
%   chancenter()

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

function E = unitsph(E, thresh)

if nargin < 2
    thresh = 1e-6;
end

[th, phi, r] = cart2sph(E(:, 1), E(:, 2), E(:, 3));

if any(abs(r - 1) > thresh)
    [E(:, 1), E(:, 2), E(:, 3)] = chancenter(E(:, 1), E(:, 2), E(:, 3), []); % Re-center
    [th, phi] = cart2sph(E(:, 1), E(:, 2), E(:, 3));
    [E(:, 1), E(:, 2), E(:, 3)] = sph2cart(th, phi, 1); % Project to unit sphere
end
