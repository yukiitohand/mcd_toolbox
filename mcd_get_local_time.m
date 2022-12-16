function [params] = mcd_get_local_time(dt,xlon,xlat,varargin)
% [params] = mcd_get_local_time(dt,xlon,xlat,varargin)
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
%   'XZ': scalar or array
%       Vertical coordinate of the requested point. Its exact definition 
%       depends on the value of input argument zkey. With the default 
%       zkey=3, z is height above surface (m).
%       (default) linspace(0,20000,1000)
%   'EXT_KEYS': cell array
%       input for mcd_set_extvarkey
%       (default) 
%            {'Ls','LTST','universal_solar_time','lmst'}
%   see mcd_query for additional optional parameters. If you directly 
%   specify 'EXTVARKEYS', then EXT_KEYS is overridden. but currently
%   extvarkeys won't be returend. not implemented yet.
%       
%% VARARGIN
xz = 0;
ext_keys = {'Ls','LTST','universal_solar_time','lmst'};
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'XZ'
                xz = varargin{i+1};
            case 'EXT_KEYS'
                ext_keys = varargin{i+1};
        end
    end
end

%% get julian date from dt (datetiem obj)
[varargin_mcd_get_julian_date] = trim_varargin_pairs(varargin,'VERBOSE');
[date_julian] = mcd_get_julian_date(dt,varargin_mcd_get_julian_date{:});

%% Do query
% set extvarkey
[extvarkey,extvarkey_idxes] = mcd_set_extvarkey(ext_keys{:});

% get information
[varargin_mcd_query] = trim_varargin_pairs(varargin,...
    {'ZKEY','HIRESKEY','DATA_SET','SCENA','EXTVARKEYS','VERBOSE'});

params = [];

[pres,rho,temp,extvar] = mcd_query(date_julian,xz,xlon,xlat,...
    'EXTVARKEYS',extvarkey,varargin_mcd_query{:});

% store information
params.xz   = xz;
params.pres = pres(1);
params.rho  = rho(1);
params.temp = temp(1);
for i=1:length(ext_keys)
    ext_key = ext_keys{i};
    params.(ext_key) = extvar(extvarkey_idxes(i));
end

end
