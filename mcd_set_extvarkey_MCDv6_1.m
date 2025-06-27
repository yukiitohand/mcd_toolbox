function [extvarkey,extvarkey_idxes] = mcd_set_extvarkey_MCDv6_1(varargin)
% [extvarkey,extvarkey_idxes] = mcd_set_extvarkey_MCDv6_1(varargin)
% Set up extvarkey for mcd_query
% INPUT Parameters
%   key1, key2, ...: string
%       key for which you want to set extvarkey=1
% OUTPUT Parameters
%   extvarkey: array, [1,100]
%       values for which key is set are ones, otherwise zeros.
%   extvarkey_idxes: array, length same as the number of keys.
%       indexes for each keys.

len_keys = length(varargin);

extvarkey = zeros(1,100);
extvarkey_idxes = zeros(1,len_keys);

for i=1:len_keys
    key = varargin{i};
    switch lower(key)
        case {'radial_dist_to_planet_center'}
            extvarkey_idxes(i) = 1;
        case {'alt_above_areoid'}
            extvarkey_idxes(i) = 2;
        case {'alt_above_surf'}
            extvarkey_idxes(i) = 3;
        case {'orographic_height'}
            extvarkey_idxes(i) = 4;
        case {'ls','solar_longitude'}
            extvarkey_idxes(i) = 9;
        case {'ltst','local_true_solar_time'}
            extvarkey_idxes(i) = 10;
        case {'universal_solar_time'}
            extvarkey_idxes(i) = 12;
        case {'lmst','local_mean_solar_time'}
            extvarkey_idxes(i) = 11;
        case {'surf_h2oice'}
            extvarkey_idxes(i) = 45;
        case {'surf_co2ice'}
            extvarkey_idxes(i) = 44;
        case {'col_h2ovapor','c_h2ovapor'}
            extvarkey_idxes(i) = 47;
        case {'vmr_h2ovapor','vol_mix_ratio_h2ovapor'}
            extvarkey_idxes(i) = 48;
        case {'col_h2oice','c_h2oice'}
            extvarkey_idxes(i) = 49;
        case {'vmr_h2oice','vol_mix_ratio_h2oice'}
            extvarkey_idxes(i) = 50;
        case {'vmr_co2','vol_mix_ratio_co2'}
            extvarkey_idxes(i) = 64;
        case {'vmr_co','vol_mix_ratio_co'}
            extvarkey_idxes(i) = 67;
        case {'c_co2','col_co2'}
            extvarkey_idxes(i) = 74;
        case {'c_co','col_co'}
            extvarkey_idxes(i) = 77;
        otherwise
            error('Key %s is not defined', key);
    end
end

extvarkey(extvarkey_idxes) = 1;

end
