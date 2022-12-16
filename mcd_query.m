function [pres,rho,temp,extvar] = mcd_query(date_julian,xz,xlon,xlat,varargin)
% [pres,rho,temp,extvar] = mcd_query(date_julian,xz,xlon,xlat,varargin)
%
% INPUT Parameters
%   dt: datetime obj
%       Earth date
%   xz: Vertical coordinate of the requested point. Its exact definition 
%       depends on the value of input argument zkey. With the default 
%       zkey=3, z is height above surface (m).
%   xlon: East Longitude (planetocentric), in degrees.
%   xlat: Latitude (planetocentric), in degrees.
%
% OUTPUT Parameters
%   pres  : scalar, atmospheric pressure (Pa)
%   rho   : scalar, atmospheric density (kg/m^3)
%   temp  : scalar, atmospheric temperature (K)
%   extvar: extra variables (array of 100)
% OPTIONAL Parameters
%   'ZKEY': integer, {1,2,3,4,5}
%       type of vertical coordinate xz
%        | 1 = radius from centre of planet (m)
%        | 2 = height above areoid (m) (MOLA zero datum)
%        | 3 = height above surface (m)
%        | 4 = pressure level (Pa)
%        |_5 = altitude above mean Mars Radius(=3396000m) (m)
%       (default) 3
%   'HIRESKEY': integer {0,1}
%       (default) 1
%   'DATA_SET': string
%       path to the data base
%       (default) '/Users/yukiitoh/data/MCD_DATA/'
%   'SCENA': integer
%       1 = Climatology ave solar
%       2 = Climatology min solar
%       3 = Climatology max solar
%       4 = dust storm tau=5 (dark dust) min solar           
%       5 = dust storm tau=5 (dark dust) ave solar           
%       6 = dust storm tau=5 (dark dust) max solar           
%       7 = warm scenario - dusty, max solar        
%       8 = cold scenario - low dust, min solar
%       24 = Mars Year 24, with associated solar EUV
%       25 = Mars Year 25, with associated solar EUV
%       26 = Mars Year 26, with associated solar EUV
%       27 = Mars Year 27, with associated solar EUV
%       28 = Mars Year 28, with associated solar EUV
%       29 = Mars Year 29, with associated solar EUV
%       30 = Mars Year 30, with associated solar EUV
%       31 = Mars Year 31, with associated solar EUV
%       32 = Mars Year 32, with associated solar EUV
%       33 = Mars Year 32, with associated solar EUV
%       (default) 1
%   'EXTVARKEYS': array [1,100]
%       (default) zeros(1,100)
%       
%   'VERBOSE': boolean
%       whether or not to print some
%       (default) 1


%% VARARGIN
zkey     = 3;
hireskey = 1;
verbose  = true;
dset     = '/project/crism/users/itohy1/data/MCD/MCD5.3/data/';
scena    = 1;
extvarkeys = zeros(1,100);
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'ZKEY'
                zkey = varargin{i+1};
            case 'HIRESKEY'
                hireskey = varargin{i+1};
            case 'DATA_SET'
                dset = varargin{i+1};
            case 'SCENA'
                scena = varargin{i+1};
                if length(scena)>50
                    error('SCENA needs to be shorter than 50 chars');
                end
            case 'EXTVARKEYS'
                extvarkeys = varargin{i+1};
                if ~isvector(extvarkeys) || length(extvarkeys) ~= 100
                    error('EXTVARKEYS needs to be 100 length vector');
                end
            case 'VERBOSE'
                verbose = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

%%
% mlb_call_mcd()
% Arguments (inputs):
%
%TIME
datekey   = 0;            % <integer> type of input date (1=Mars date, 0=Earth date) 
localtime = 0;            % <real> local time (in martian hr) at lon. xlon 
                          % if datekey=1, then should be set to 0.
xdate     = date_julian;  % <double precision> date 
                          % (IF datekey = 1 : Value of Ls [deg.])
                          % (IF datekey = 0 : julian date)
% LOCATION
% zkey      = 3;            % <integer> type of vertical coordinate xz (3 = above surface [m])
% xz        = crd(3);
% xz        = 100.0;      % <real> vertical coordinate (m or Pa, depends on zkey)
% set hires flag
% hireskey  = 1;            % <integer> (1 = switch to high res. topography)
% xlat      = crd(1);
% xlon      = crd(2);
% xlat      = 5.0;        % <real> latitude (degrees north)
% xlon      = 6.0;        % <real> longitude (degrees east)

% DUST scenario
% dset      = 'MCD_DATA/'; % <character*50> data set
% scena     = 1;           % <integer> scenario (1 = Climatology ave solar)
perturkey = 1;           % <integer>  perturbation type (1= none)
seedin    = 1;           % <real> 
gwlength  = 0.0;         % <real>  for small scale (ie: gravity wave) perturbations;
% extvarkeys(1)=0;         % <integer> array output type (extvar(i) = 0 : don't compute)

% for i=2:100
%  extvarkeys(i)=extvarkeys(1);
% end % i
%
% Output variables must be mentioned before the function call in order to 
% reserve memory 
% Arguments (outputs):
%
meanvar = zeros(1,5);  % <real> mean unperturbed values (array of 5)
extvar = zeros(1,100); % <real>  extra variables (array of 100)

pres    = zeros(1,2); % <real> atmospheric pressure (Pa)
rho    = zeros(1,2); % <real> atmospheric density (kg/m^3)
temp    = zeros(1,2); % <real> atmospheric temperature (K)
u       = zeros(1,2); % <real> zonal wind component (East-West)
v       = zeros(1,2);  % <real> meridional wind component (North-South)
seedout = zeros(1,2); % <real> current value of the seed of the random number generator
ier     = zeros(1,2);  % <integer> error flag (0 = OK)

mlb_call_mcd(zkey,xz,xlon,xlat,hireskey,datekey,xdate,localtime,dset,scena,...
    perturkey,seedin,gwlength,extvarkeys,pres,rho,temp,u,v,meanvar,extvar,seedout,ier);

if verbose
    if ier(1)==0
        fprintf('p   = %g Pa\n', pres(1));
        fprintf('rho = %g kg/m3\n',rho(1));
        fprintf('T   = %g K \n',temp(1));
        fprintf('Zonal wind      = %g m/s\n',u(1));
        fprintf('Meridional wind = %g m/s\n',v(1));
        fprintf('\n');
        for i=1:5
            fprintf('meanvar(%d)= %g\n',i,meanvar(i));
        end
        fprintf('\n'); 
        for i=1:7
            fprintf('extvar(%d) = %g\n',i,extvar(i));
        end
        for i = 8:80
            if extvarkeys(i)~=0 
                fprintf('extvar(%d) = %g\n',i,extvar(i));
            end
        end
    else
        fprintf(2, 'CALL_MCD ERROR !!\n');
        fprintf(2, 'returned error code: %d\n',ier(1));
    end
else
    if ier(1)~=0
        fprintf(2, 'CALL_MCD ERROR !!\n');
        fprintf(2, 'returned error code: %d\n',ier(1));
    end
end

ires = 0;

end

