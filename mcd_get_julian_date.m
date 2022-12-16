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
%   'VERBOSE': boolean
%       whether or not to print some
%       (default) 1
%% INPUT CHECK
if ~isdatetime(dt)
    error('Input dt is not a datetime obj.');
end

%% VARARGIN
verbose = true;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
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

% Output variables must be mentioned before the function call in order to 
% reserve memory 
% Arguments (outputs):
ierr = zeros(1,2);           % <integer> error flag (0 = OK)
date_julian = zeros(1,2);    % <real> julian date            

mlb_julian(month,day,year,hour,minute,second,ierr,date_julian);

if verbose
    if ierr(1)==0
     fprintf('month  = %g\n', month );
     fprintf('day    = %g\n', day   );
     fprintf('year   = %g\n', year  );
     fprintf('hour   = %g\n', hour  );
     fprintf('minute = %g\n', minute);
     fprintf('second = %g\n', second);
     fprintf('This is julian date = %g\n\n', date_julian(1));
    else
     fprintf('julian ERROR !!\n');
     fprintf('returned error code: %d\n\n',ierr(1));
    end
end

date_julian = date_julian(1);

end