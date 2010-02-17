% eegplugin_sphspline() - EEGLAB plugin providing basic spherical spline
%                         interpolation functions
%
% Usage:
%   >> eegplugin_sphspline(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Author: Andreas Widmann, University of Leipzig, Germany, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

% $Log: eegplugin_sphspline.m,v $

function vers = eegplugin_sphspline(fig, trystrs, catchstrs)

    vers = 'sphspline0.2beta2';
    if nargin < 3
        error('eegplugin_sphspline requires 3 arguments');
    end

    % add folder to path
    if ~exist('sserp.m')
        p = which('eegplugin_sphspline');
        p = p(1:findstr(p,'eegplugin_sphspline.m')-1);
        addpath([p vers]);
    end

    % menu callbacks
    comrepchan      = [trystrs.no_check '[EEG LASTCOM] = pop_repchan(EEG);'      catchstrs.new_and_hist];
    comscdtransform = [trystrs.no_check '[EEG LASTCOM] = pop_scdtransform(EEG);' catchstrs.new_and_hist];
    comaddchanica   = [trystrs.no_check '[EEG LASTCOM] = pop_addchanica(EEG);'   catchstrs.new_and_hist];
    complotsserpmap = [trystrs.no_check '[LASTCOM] = pop_plotsserpmap(EEG);'     catchstrs.add_to_hist];

    % create submenus and items
    menu = findobj(fig, 'tag', 'tools');
    submenu = uimenu(menu, 'Label', 'Spherical spline interpolation');
    uimenu(submenu, 'Label', 'Replace bad channel', 'CallBack', comrepchan);
    uimenu(submenu, 'Label', 'SCD transform the data', 'CallBack', comscdtransform);
    uimenu(submenu, 'Label', 'Add channel to ICA', 'CallBack', comaddchanica);
    menu  = findobj(fig, 'tag', 'plot');
    submenu = uimenu(menu, 'Label', 'Spherical spline interpolation');
    uimenu(submenu, 'Label', '2D map series', 'CallBack', complotsserpmap);
