% pop_plotsserpmap() - Plot 2D map series using spherical spline
%                      interpolation
%
% Usage:
%   >> [EEG, com] = pop_plotsserpmap(EEG); % pop-up window mode
%   >> [EEG, com] = pop_plotsserpmap(EEG, 'parameter1', value1, ...
%                                         'parameter2', value2, ...
%                                         'parametern', valuen);
%
% Inputs:
%   EEG           - EEGLAB EEG structure
%
% Optional inputs:
%   'plot'        - string plot type. 'erp' (plot channel ERP maps) or
%                   'comp' (plot component maps) {default 'erp'}
%   'type'        - string type of interpolation. 'sp' (scalp
%                   potential), 'scd' (scalp current density), or 'lap'
%                   (surface Laplacian) {default 'sp'}
%   'items'       - scalar, vector, or 2 by n matrix items to plot. Time
%                   points (s) (window start and end times in rows) or
%                   component numbers
%   'lambda'      - scalar smoothing factor (commonly used values are
%                   1e-7 for sp and 1e-5 for scd) {default 0}
%   'proj'        - string type of projection. 'stereographic',
%                   'orthographic', 'gnomonic', 'equidistant', or
%                   'equiareal' {default 'equiareal'}
%   'rowsCols'    - vector of integers [rows columns]
%   'colormap'    - string colormap {default 'jet'}
%   'maplimits'   - string or 2 element vector maplimits. 'absmax'
%                   (zero-symmetric color mapping), 'maxmin' (color
%                   mapping scaled to global extrema), or [min max]
%                   {default 'absmax'}
%
% Outputs:
%   com           - history string
%
% Note:
%   Channels without location coordinates are removed. Channel
%   coordinates should be located on the surface of a (unit) sphere.
%
% Author: Andreas Widmann, University of Leipzig, 2006
%
% See also:
%   sserpgfcn(), sserpweights(), sserp(), unitsph(), plotsserpmap(),
%   plotsserpmapproj()

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

function [com] = pop_plotsserpmap(EEG, varargin)

com = '';

% Arguments
if nargin < 1
    help pop_plotsserpmap
    return
end
if isempty(EEG.data)
    error('Cannot process empty dataset.')
end

% Defaults
COLOR_EEGLAB_BACKGROUND = [.93 .96 1];

if nargin < 2
    drawnow;
    [Plot(1:2).name] = deal('ERP maps', 'Component maps');
    [Plot(1:2).arg] = deal('erp', 'comp');
    [Type(1:3).name] = deal('Scalp potential', 'Scalp current density', 'Surface Laplacian');
    [Type(1:3).arg] = deal('sp', 'scd', 'lap');
    [Proj(1:5).name] = deal('Lambert azimuthal equiareal', 'Azimuthal equidistant', 'Orthographic', 'Stereographic', 'Gnomonic');
    [Proj(1:5).arg] = deal('equiareal', 'equidistant', 'orthographic', 'stereographic', 'gnomonic');
    uigeom = {[1 1] [1 1] [1 1] [1 1] [1 1] [1 1]};
    uilist = {{'style' 'text' 'string' 'Plot:'} ...
              {'style' 'popupmenu' 'string' {Plot.name} 'value' 1} ...
              {'style' 'text' 'string' 'Type:'} ...
              {'style' 'popupmenu' 'string' {Type.name} 'value' 1} ...
              {'style' 'text' 'string' 'Window start times (ms)/ components:'} ...
              {'style' 'edit' 'string' ''} ...
              {'style' 'text' 'string' 'Window end times (ms):'} ...
              {'style' 'edit' 'string' ''} ...
              {'style' 'text' 'string' 'Smoothing {default 0}:'} ...
              {'style' 'edit' 'string' ''} ...
              {'style' 'text' 'string' 'Projection:'} ...
              {'style' 'popupmenu' 'string' {Proj.name} 'value' 1}};
    result = inputgui(uigeom, uilist, 'pophelp(''pop_plotsserpmap'')', 'Plot 2D map series -- pop_plotsserpmap()');
    if isempty(result), return, end

    Arg = {'plot' Plot(result{1}).arg 'type' Type(result{2}).arg 'proj' Proj(result{6}).arg};
    Arg = [Arg {'items'} {str2num(result{3})}];
    Arg{end} = [Arg{end}; str2num(result{4})];
    Arg = [Arg {'lambda'} {str2double(result{5})}];
else
    Arg = varargin;
end

% Arguments and defaults
Arg = struct(Arg{:});
if ~isfield(Arg, 'items') || isempty(Arg.items)
    error('Not enough input arguments.')
end
if ~isfield(Arg, 'plot') || isempty(Arg.plot)
    Arg.plot = 'erp';
end
if ~isfield(Arg, 'type') || isempty(Arg.type)
    Arg.type = 'sp';
end
if ~isfield(Arg, 'lambda') || isempty(Arg.lambda)
    Arg.lambda = 0;
end
if ~isfield(Arg, 'proj') || isempty(Arg.proj)
    Arg.proj = 'equiareal';
end
if ~isfield(Arg, 'rowsCols') || isempty(Arg.rowsCols)
    Arg.rowsCols(2) = ceil(sqrt(size(Arg.items, 2) + double(size(Arg.items, 2) > 1)));
    Arg.rowsCols(1) = ceil((size(Arg.items, 2) + double(size(Arg.items, 2) > 1)) / Arg.rowsCols(2));
end
if ~isfield(Arg, 'colormap') || isempty(Arg.colormap)
    Arg.colormap = 'jet';
end
if ~isfield(Arg, 'maplimits') || isempty(Arg.maplimits)
    Arg.maplimits = 'absmax';
end
Units.erp.sp = 'µV';
Units.erp.scd = 'mA/m^3';
Units.erp.lap = 'mV/m^2';
Units.comp.sp = '';
Units.comp.scd = '';
Units.comp.lap = '';

% ms to sample
if strcmpi(Arg.plot, 'erp')
    pntArray = round((Arg.items / 1000 - EEG.xmin) * EEG.srate + 1);
    Arg.items = ((pntArray  - 1) / EEG.srate + EEG.xmin) * 1000;
    if size(pntArray, 1) == 1
        pntArray(2, :) = pntArray(1, :);
    end
end

% Channels with location coordinates
chanArray = ~sum(cellfun('isempty', {EEG.chanlocs.X; EEG.chanlocs.Y; EEG.chanlocs.Z}));

% Unit sphere
E = unitsph([EEG.chanlocs(chanArray).X; EEG.chanlocs(chanArray).Y; EEG.chanlocs(chanArray).Z]');

% Collect arguments for interpolation
[F, ePrime, x, convexHull, headRadius] = plotsserpmapproj(E, Arg.proj);
isErp = ~sum(isnan(F), 2); % NaNs slow down interpolation by factor 10 to 100
[Ginv, g] = sserpgfcn(E, F(isErp, :), Arg.type, Arg.lambda);
erpData = NaN(length(x), length(x));

% Open figure
hfig = figure;
set(hfig, 'color', COLOR_EEGLAB_BACKGROUND, 'Colormap', colormap(Arg.colormap))

zMin = []; zMax = [];

for iAxis = 1 : size(Arg.items, 2);

    if strcmpi(Arg.plot, 'erp')
        % Average data
        data = mean(mean(EEG.data(chanArray, pntArray(1, iAxis) : pntArray(2, iAxis), :), 2), 3);
    else
        data = EEG.icawinv(chanArray, Arg.items(iAxis));
    end

    % Interpolate
    C = sserpweights(data, Ginv);
    erpData(isErp) = sserp(C, g, Arg.type);

    % Find extrema
    zMin = [zMin min(min(erpData))];
    zMax = [zMax max(max(erpData))];
    if strcmpi(Arg.plot, 'comp')
        zAbsMax = max(abs([zMin(end) zMax(end)]));
        erpData = erpData / zAbsMax;
        zMin(end) = zMin(end) / zAbsMax;
        zMax(end) = zMax(end) / zAbsMax;
    end

    % Plot map
    subplot(Arg.rowsCols(1), Arg.rowsCols(2), iAxis)
    h(iAxis) = plotsserpmap(erpData, ePrime, x, convexHull, headRadius);

end

% Scale maps to maplimits
if strcmpi(Arg.maplimits, 'absmax')
    Arg.maplimits = [-max(abs([zMin zMax])) max(abs([zMin zMax]))];
elseif strcmpi(Arg.maplimits, 'maxmin')
    Arg.maplimits = [min(zMin) max(zMax)];
end
set([h.axis], 'Clim', Arg.maplimits);

% Colorbar
if size(Arg.items, 2) > 1
    subplot(Arg.rowsCols(1), Arg.rowsCols(2), Arg.rowsCols(1) * Arg.rowsCols(2), 'Visible', 'off');
    h(end + 1).axis = colorbar('v6');
    set(get(h(end).axis, 'Children'), 'YData', [Arg.maplimits(1) + diff(Arg.maplimits) / 128 Arg.maplimits(2) - diff(Arg.maplimits) / 128])
    set(h(end).axis, 'YLim', Arg.maplimits)
    drawnow
else
    h(end + 1).axis = colorbar('v6');
end
set(h(end).axis, 'YTickMode', 'manual', 'NextPlot', 'Add')

% Adjust contour line levels to colorbar ticks
levelList = get(h(end).axis, 'YTick');
if length(levelList(levelList > min(zMin) & levelList < max(zMax))) < 4
    tickStep = diff(levelList(1:2)) / 4;
    levelList = levelList(1) - 3 * tickStep : tickStep : levelList(end) + 3 * tickStep;
elseif length(levelList(levelList > min(zMin) & levelList < max(zMax))) < 6
    tickStep = diff(levelList(1:2)) / 2;
    levelList = levelList(1) - tickStep : tickStep : levelList(end) + tickStep;
end
levelList = levelList(levelList > min(zMin) & levelList < max(zMax));

% Set contour lines
set([h(1:end - 1).contour], 'LevelList', levelList);

% Add contour lines to colorbar
axes(h(end).axis);
[foo, h(end).contour] = contour(get(h(end).axis, 'XLim'), [Arg.maplimits], [Arg.maplimits' Arg.maplimits'], levelList, 'LineColor', [0 0 0]);

% Captions
if strcmpi(Arg.plot, 'erp')
    if size(Arg.items, 1) == 1
        captions = cellstr(num2str(Arg.items', '%.0f ms'));
    else
        captions = cellstr(num2str(Arg.items', '%.0f - %.0f ms'));
    end
else
    captions = cellstr(num2str(Arg.items', '%d'));
end
title = get([h(1:end - 1).axis], 'Title');
if ~iscell(title), title = {title}; end
set([title{:}], {'String'}, captions, 'Visible', 'on')

% Units
xLabel = get([h(end).axis], 'XLabel');
set(xLabel, 'String', Units.(Arg.plot).(Arg.type), 'Visible', 'on')

% History string
com = ['pop_plotsserpmap(' inputname(1) ', ' arg2str(Arg) ');'];
