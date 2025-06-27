function [params] = mcd_get_vertical_profile(dt,xlon,xlat,varargin)
% [params] = mcd_get_vertical_profile(dt,xlon,xlat,ext_keys,varargin)
% 
% INPUT Parameters
%   dt: datetime obj, earth time
%   xlon: East Longitude (planetocentric), in degrees.
%   xlat: Latitude (planetocentric), in degrees.
% OUTPUT Parameters
%   params: struct having fields,
%       pres  : scalar, atmospheric pressure (Pa)
%       rho   : scalar, atmospheric density (kg/m^3)
%       temp  : scalar, atmospheric temperature (K)
%       ext_keys{i}: scalar, any extra parameters
% OPTIONAL Parmeters
%   'MCD_VER': str ['6_1', '5_3']
%       Version of MCD
%       (default) '6_1'
%   'XZ': scalar or array
%       Vertical coordinate of the requested point. Its exact definition 
%       depends on the value of input argument zkey. With the default 
%       zkey=3, z is height above surface (m).
%       (default) linspace(0,20000,1000)
%   'EXT_KEYS': cell array
%       input for mcd_set_extvarkey
%       (default) 
%            {'radial_dist_to_planet_center','alt_above_areoid',...
%             'alt_above_surf','orographic_height',...
%             'Ls','LTST','universal_solar_time',...
%             'surf_h2oice','surf_co2ice','vmr_h2oice','col_h2oice',...
%             'vol_mix_ratio_co2','vol_mix_ratio_co','vol_mix_ratio_h2ovapor',...
%             'col_co2','col_co','col_h2ovapor'}
%   see mcd_query for additional optional parameters. If you directly 
%   specify 'EXTVARKEYS', then EXT_KEYS is overridden. but currently
%   extvarkeys won't be returend. not implemented yet.
%       
%% VARARGIN
xz = linspace(0,50000,1000);
ext_keys = {'radial_dist_to_planet_center','alt_above_areoid',...
            'alt_above_surf','orographic_height',...
            'Ls','LTST','universal_solar_time',...
            'surf_h2oice','surf_co2ice','vmr_h2oice','col_h2oice',...
            'vol_mix_ratio_co2','vol_mix_ratio_co','vol_mix_ratio_h2ovapor',...
            'col_co2','col_co','col_h2ovapor'};
mcd_ver = '6_1';
scena = 1;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'XZ'
                xz = varargin{i+1};
            case 'SCENA'
                scena = varargin{i+1};
            case 'MCD_VER'
                mcd_ver = varargin{i+1};
            case 'EXT_KEYS'
                ext_keys = varargin{i+1};
        end
    end
end

%% get julian date from dt (datetiem obj)
[varargin_mcd_get_julian_date] = trim_varargin_pairs(varargin, ...
    {'VERBOSE', 'MCD_VER'});
[date_julian] = mcd_get_julian_date(dt,varargin_mcd_get_julian_date{:});

%% Do query
% set extvarkey
switch mcd_ver
    case '6_1'
        [extvarkey,extvarkey_idxes] = ...
            mcd_set_extvarkey_MCDv6_1(ext_keys{:});
    case '5_3'
        [extvarkey,extvarkey_idxes] = ...
            mcd_set_extvarkey_MCDv5_3(ext_keys{:});
    otherwise
        error('MCD_VER must be either ''6_1'' or ''5_3''.');
end

% get information
[varargin_mcd_query] = trim_varargin_pairs(varargin,...
    {'ZKEY', 'HIRESKEY', 'DATA_SET', 'SCENA', 'EXTVARKEYS', 'VERBOSE', ...
    'MCD_VER'});

params = [];
for n=1:length(xz)
    xzi = xz(n);
    [pres,rho,temp,extvar] = mcd_query(date_julian,xzi,xlon,xlat,...
        'EXTVARKEYS',extvarkey,varargin_mcd_query{:});

    % store information
    params(n).xz   = xzi;
    params(n).pres = pres(1);
    params(n).rho  = rho(1);
    params(n).temp = temp(1);
    for i=1:length(ext_keys)
        ext_key = ext_keys{i};
        params(n).(ext_key) = extvar(extvarkey_idxes(i));
    end
end

end
