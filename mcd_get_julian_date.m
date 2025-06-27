function [date_julian] = mcd_get_julian_date(dt,varargin)
% [date_julian] = mcd_get_julian_date(dt,varargin)
% INPUT Parameters
%   dt: datetime obj
%       Earth date
%   crd: array, [lat,lon,z]
%       lat - Latitude (planetocentric), in degrees.
%       lon - East Longitude (planetocentric), in degrees.
%       z   - Vertical coordinate of the requested point. Its exact 
%             definition depends on the value of input argument zkey.
%             With the default zkey=3, z is height above surface (m).
% OUTPUT Parameters
%   date_julian: julian date.
%  
% OPTIONAL Parameters
%   'MCD_VER': str ['6_1', '5_3']
%       Version of MCD
%       (default) '6_1'
%   'VERBOSE': boolean
%       whether or not to print some
%       (default) 1
%% INPUT CHECK
if ~isdatetime(dt)
    error('Input dt is not a datetime obj.');
end

%% VARARGIN
mcd_ver = '6_1';
verbose = true;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'MCD_VER'
                mcd_ver = varargin{i+1};
            case 'VERBOSE'
                verbose = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

%% mlb_julian
% Arguments (inputs): 
month     = dt.Month;      % <integer> month
day       = dt.Day;        % <integer> day
year      = dt.Year;       % <integer> year
hour      = dt.Hour;       % <integer> hour
minute    = dt.Minute;     % <integer> minute
second    = dt.Second;     % <integer> second

switch mcd_ver
    case '6_1'
        [date_julian, ierr] = mlb_julian_MCDv6_1( ...
            month, day, year, hour, minute, second);
    case '5_3'
        [date_julian, ierr] = mlb_julian_MCDv5_3( ...
            month, day, year, hour, minute, second);
    otherwise
        error('Supported MCD_VER are ''6_1'' and ''5_3''.');
end

% Two outputs are returned.
% date: <real> julian date
% ierr: <integer> error flag (0 = OK)


if ierr(1)==0
    if verbose
        fprintf('month  = %g\n', month );
        fprintf('day    = %g\n', day   );
        fprintf('year   = %g\n', year  );
        fprintf('hour   = %g\n', hour  );
        fprintf('minute = %g\n', minute);
        fprintf('second = %g\n', second);
        fprintf('This is julian date = %g\n\n', date_julian(1));
    end
else
    fprintf('julian ERROR !!\n');
    fprintf('returned error code: %d\n\n',ierr(1));
end

end