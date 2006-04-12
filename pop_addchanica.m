% pop_addchanica() - Add channel to icawinv and back project to data
%
% Usage:
%   >> [EEG, com] = pop_addchanica(EEG); % pop-up window mode
%   >> [EEG, com] = pop_addchanica(EEG, 'parameter1', value1, ...
%                                       'parameter2', value2, ...
%                                       'parametern', valuen);
%
% Inputs:
%   EEG       - EEGLAB EEG structure
%   'nums'    - scalar or vector new channel numbers
%   'labels'  - cell array of strings new channel labels
%   'locs'    - string file with channel location coordinates
%
% Optional inputs:
%   'nTerms'  - scalar int > 0 number of terms {default 50}
%   'm'       - scalar int > 1 m {default 4}
%
% Outputs:
%   EEG       - EEGLAB EEG structure
%   com       - history string
%
% Note:
%   Channel coordinates should be located on the surface of a (unit)
%   sphere.
%
% Author: Andreas Widmann, University of Leipzig, 2006
%
% See also:
%   sserp(), sserpgfcn(), sserpweights(), unitsph(), pop_chanedit()

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

function [EEG, com] = pop_addchanica(EEG, varargin)

com = '';

if nargin < 1
    help pop_addchanica;
    return
end
if isempty(EEG.icawinv)
    error('Cannot process empty dataset.');
end

% GUI
if nargin < 2
    drawnow;
    uigeom = {[1 1 0.5] [1 1 0.5] [1 1 0.5] [1 0.5 1] [1 0.5 1]};
    uilist = {{'style' 'text' 'string' 'New channel number(s):'} ...
              {'style' 'edit' 'string' ''} {} ...
              {'style' 'text' 'string' 'New channel label(s):'} ...
              {'style' 'edit' 'string' ''} {} ...
              {'style' 'text' 'string' 'File with channel location coordinates:'} ...
              {'style' 'edit' 'string' '' 'tag' 'locsedit'} ...
              {'style' 'pushbutton' 'string' 'Browse' 'callback' '[filename pathname] = uigetfile2(''*.*'', ''Load a channel location file''); set(findobj(gcbf, ''tag'', ''locsedit''), ''string'', fullfile(pathname, filename));'} ...
              {'style' 'text' 'string' 'Number of terms {default 50}:'} ...
              {'style' 'edit' 'string' ''} {} ...
              {'style' 'text' 'string' 'm {default 4}:'} ...
              {'style' 'edit' 'string' ''} {}};
    result = inputgui(uigeom, uilist, 'pophelp(''addicachan'')', 'Add channel to ICA -- pop_addchanica()');
    if isempty(result), return; end

    Arg = {};
    if ~isempty(result{1})
        Arg = [Arg {'nums'} {str2num(result{1})}];
    end
    if ~isempty(result{2})
        Arg = [Arg {'labels'} {eval(['{' result{2} '}'])}];
    end
    if ~isempty(result{3})
        Arg = [Arg {'locs'} {result{3}}];
    end
    if ~isempty(result{4})
        Arg = [Arg {'nTerms'} {str2double(result{2})}];
    end
    if ~isempty(result{5})
        Arg = [Arg {'m'} {str2double(result{3})}];
    end
else
    Arg = varargin;
end

% Check parameters
Arg = cell2struct(Arg(2:2:end), Arg(1:2:end - 1), 2);
if ~all(isfield(Arg, {'nums' 'labels' 'locs'})) || isempty(Arg.nums) || isempty(Arg.labels) || isempty(Arg.locs)
    error('Not enough input arguments.');
end
if ~isfield(Arg, 'nTerms') || isempty(Arg.nTerms)
    Arg.nTerms = 50;
end
if ~isfield(Arg, 'm') || isempty(Arg.m)
    Arg.m = 4;
end

% Insert channels
chanArray = ~ismember(1 : EEG.nbchan + length(Arg.nums), Arg.nums);
EEG.nbchan = length(chanArray);
EEG.icawinv(chanArray, :) = EEG.icawinv;
EEG.data(chanArray, :, :) = EEG.data;

% Lookup channel location coordinates
Chanlocs(chanArray) = EEG.chanlocs;
[Chanlocs(Arg.nums).labels] = deal(Arg.labels{:});
EEG.chanlocs = pop_chanedit(Chanlocs, 'lookup', Arg.locs);

% Find channels with location coordinates
chanArray = chanArray & ~sum(cellfun('isempty', {EEG.chanlocs.X; EEG.chanlocs.Y; EEG.chanlocs.Z}));

% Channel location coordinate matrices E and F and unit sphere
E = unitsph([EEG.chanlocs(chanArray).X; EEG.chanlocs(chanArray).Y; EEG.chanlocs(chanArray).Z]');
F = unitsph([EEG.chanlocs(Arg.nums).X; EEG.chanlocs(Arg.nums).Y; EEG.chanlocs(Arg.nums).Z]');
if size(F, 1) ~= length(Arg.nums)
    error('No channel location coordinates for channel(s) to interpolate.');
end

% G matrix
[Ginv, g] = sserpgfcn(E, F, 'sp', 0, Arg.nTerms, Arg.m);

% Interpolate
C = sserpweights(EEG.icawinv(chanArray, :), Ginv);
EEG.icawinv(Arg.nums, :) = sserp(C, g, 'sp');
EEG.icasphere = pinv(EEG.icaweights) * pinv(EEG.icawinv);

% Back project new channel to EEG.data
EEG.data(Arg.nums, :, :) = reshape(EEG.icawinv(Arg.nums, :) * EEG.icaact(:, :), [length(Arg.nums) EEG.pnts EEG.trials]);

% History string
com = [inputname(1) ' = pop_addchanica(' inputname(1) ', ' arg2str(Arg) ');'];
