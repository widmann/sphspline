% pop_repchan() - Replace bad channel(s) using spherical spline
%                 interpolation
%
% Usage:
%   >> [EEG, com] = pop_repchan(EEG); % pop-up window mode
%   >> [EEG, com] = pop_repchan(EEG, 'parameter1', value1, ...
%                                    'parameter2', value2, ...
%                                    'parametern', valuen);
%
% Inputs:
%   EEG       - EEGLAB EEG structure
%   'chans'   - scalar or vector channel(s) to replace
%
% Optional inputs:
%   'nTerms'  - scalar int > 0 number of terms {default 50}
%   'm'       - scalar int > 1 m {default 4}
%   'lookup'  - string file with channel location coordinates
%   'nFrames' - scalar int > 0 number of frames to interpolate
%               vectorized (trade-off between speed and memory usage)
%               {default 1000}
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

function [EEG, com] = pop_repchan(EEG, varargin)

com = '';

if nargin < 1
    help pop_repchan;
    return
end
if isempty(EEG.data)
    error('Cannot process empty dataset.');
end

% GUI
if nargin < 2
    drawnow;
    uigeom = {[1 1 0.5] [1 0.5 1] [1 0.5 1] [1 1 0.5]};
    uilist = {{'style' 'text' 'string' 'Channels to replace:'} ...
              {'style' 'edit' 'string' '' 'tag' 'chanEdit'} ...
              {'style' 'pushbutton' 'string' 'Select' 'callback' 'set(findobj(gcbf, ''tag'', ''chanEdit''), ''string'', num2str(pop_chansel({EEG.chanlocs.labels}, ''withindex'', ''on'', ''selectionmode'', ''multiple'')));'} ...
              {'style' 'text' 'string' 'Number of terms {default 50}:'} ...
              {'style' 'edit' 'string' ''} {} ...
              {'style' 'text' 'string' 'm {default 4}:'} ...
              {'style' 'edit' 'string' ''} {} ...
              {'style' 'text' 'string' 'Look up channel locations from file:'} ...
              {'style' 'edit' 'string' '' 'tag' 'lookupEdit'} ...
              {'style' 'pushbutton' 'string' 'Browse' 'callback' '[filename pathname] = uigetfile2(''*.*'', ''Load a channel location file''); set(findobj(gcbf, ''tag'', ''lookupEdit''), ''string'', fullfile(pathname, filename));'}};
    result = inputgui(uigeom, uilist, 'pophelp(''pop_repchan'')', 'Replace bad channel(s) -- pop_repchan()');
    if isempty(result), return; end

    Arg = {};
    if ~isempty(result{1})
        Arg = [Arg {'chans'} {str2num(result{1})}];
    end
    if ~isempty(result{2})
        Arg = [Arg {'nTerms'} {str2double(result{2})}];
    end
    if ~isempty(result{3})
        Arg = [Arg {'m'} {str2double(result{3})}];
    end
    if ~isempty(result{4})
        Arg = [Arg {'lookup'} {result{4}}];
    end
else
    Arg = varargin;
end

% Check parameters
Arg = struct(Arg{:});
if ~isfield(Arg, 'chans') || isempty(Arg.chans)
    error('Not enough input arguments.');
end
if ~isfield(Arg, 'nTerms') || isempty(Arg.nTerms)
    Arg.nTerms = 50;
end
if ~isfield(Arg, 'm') || isempty(Arg.m)
    Arg.m = 4;
end
if ~isfield(Arg, 'nFrames') || isempty(Arg.nFrames)
    Arg.nFrames = 1000;
end

% Look up channel location coordinates
if isfield(Arg, 'lookup') && ~isempty(Arg.lookup)
    Chanlocs = pop_chanedit(EEG.chanlocs, 'lookup', Arg.lookup);
else
    Chanlocs = EEG.chanlocs;
end

% Find channels with location coordinates
chanArray = find(~sum(cellfun('isempty', {Chanlocs.X; Chanlocs.Y; Chanlocs.Z})));

% Unit sphere
E = unitsph([Chanlocs(chanArray).X; Chanlocs(chanArray).Y; Chanlocs(chanArray).Z]');

% Channel location coordinate matrices E and F
F = E(ismember(chanArray, Arg.chans), :);
if size(F, 1) ~= length(Arg.chans)
    error('No channel location coordinates for channel(s) to interpolate.');
end
E = E(~ismember(chanArray, Arg.chans), :);
chanArray = setdiff(chanArray, Arg.chans);

% G matrix
[Ginv, g] = sserpgfcn(E, F, 'sp', 0, Arg.nTerms, Arg.m);

% Reshape EEG.data
if EEG.trials > 1
    EEG.data = reshape(EEG.data, EEG.nbchan, []);
end

% Initialize progress bar
nProgBarSteps = 20;
progBarArray = ceil(linspace(size(EEG.data, 2) / nProgBarSteps, size(EEG.data, 2), nProgBarSteps));
progBarHandle = waitbar(0, '0% done', 'Name', 'Replace bad channel(s) -- pop_repchan()');
tic

blockArray = [1 : Arg.nFrames : size(EEG.data, 2) size(EEG.data, 2) + 1];
for iBlock = 1 : length(blockArray) - 1

    % C matrix
    C = sserpweights(EEG.data(chanArray, blockArray(iBlock) : blockArray(iBlock + 1) - 1), Ginv);

    % Interpolation
    EEG.data(Arg.chans, blockArray(iBlock) : blockArray(iBlock + 1) - 1) = sserp(C, g, 'sp');

    % Update progress bar
    if blockArray(iBlock + 1) - 1 >= progBarArray(1)
        progBarArray(1) = [];
        p = (nProgBarSteps - length(progBarArray)) / nProgBarSteps;
        waitbar(p, progBarHandle, [num2str(p * 100) '% done, ' num2str(ceil((1 - p) / p * toc)) ' s left']);
    end

end

% Deinitialize progress bar
if exist('progBarHandle', 'var')
    close(progBarHandle)
end

% Reshape EEG.data
if EEG.trials > 1
    EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
end

% History string
com = sprintf('%s = pop_repchan(%s', inputname(1), inputname(1));
for c = fieldnames(Arg)'
    if ischar(Arg.(c{:}))
        com = [com sprintf(', ''%s'', ''%s''', c{:}, Arg.(c{:}))];
    else
        com = [com sprintf(', ''%s'', %s', c{:}, mat2str(Arg.(c{:})))];
    end
end
com = [com ');'];
