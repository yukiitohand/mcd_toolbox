function [out] = mcd_crism_create_profile(obs_id, varargin)
% OPTIONAL Parmeters
%   'MCD_VER': str ['6_1', '5_3']
%       Version of MCD
%       (default) '6_1'

save_file = true;
mcd_ver = '6_1';
scena = 1;
save_dir = '/Volumes/LaCie5TB/out/mcd';
verbose = false;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            % ## I/O OPTIONS #---------------------------------------------
            case 'SAVE_DIR'
                save_dir = varargin{i+1};
            case 'SAVE'
                save_file = varargin{i+1};
            
            % ## MCD OPTIONS #---------------------------------------------
            case 'SCENA'
                scena = varargin{i+1};
            case 'MCD_VER'
                mcd_ver = varargin{i+1};
            
            % ## OTHER OPTIONS #-------------------------------------------
            case 'VERBOSE'
                verbose = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

% obs_id = 'C76F';
crism_obs = CRISMObservation(obs_id,'sensor_id','L');
% crism_obs = CRISMObservation('94B5','SENSOR_ID','L','OBS_COUNTER_SCENE','0[02]{1}','OBS_COUNTER_df','0[13]{1}');

switch upper(crism_obs.info.obs_classType)
    case {'FFC'}
        basenameIF = crism_obs.info.basenameIF{1};
        basenameRA = crism_obs.info.basenameRA{1};
        basenameDDR = crism_obs.info.basenameDDR{1};
    case {'FRT','HRL','FRS','HRS','ATO'}
        basenameIF = crism_obs.info.basenameIF;
        basenameDDR = crism_obs.info.basenameDDR;
        basenameRA = crism_obs.info.basenameRA;
    otherwise
        error('Not implemented for %s',upper(crism_obs.info.obs_classType));
end

TRRIFdata = CRISMdata(basenameIF,'');
sclk_img = (TRRIFdata.get_sclk_start()+TRRIFdata.get_sclk_stop())/2;
%
dt = datetime('1980-01-01 00:00:00.000') + duration(0,0,sclk_img);
DEdata = CRISMDDRdata(basenameDDR,'');
DEdata.readimg();
%% Measure time accurately
% naif_archive_init;
% SPICEMetaKrnlsObj = MRO_CRISM_SPICE_META_KERNEL(DEdata);
% SPICEMetaKrnlsObj.set_kernel_sclk(SPICEMetaKrnlsObj.diropt.sclk,'FileName',SPICEMetaKrnlsObj.src.sclk,'dwld',0);
% SPICEMetaKrnlsObj.set_kernel_lsk(SPICEMetaKrnlsObj.diropt.lsk,'FileName',SPICEMetaKrnlsObj.src.lsk,'dwld',0);
% SPICEMetaKrnlsObj.set_kernel_spk_de_default();
% SPICEMetaKrnlsObj.sclk.furnsh();
% SPICEMetaKrnlsObj.lsk.furnsh();
% SPICEMetaKrnlsObj.spk_de.furnsh();
% SC  = -74999; 
% sclk_str = spice_sclk_num2str(sclk_img);
% et  = cspice_scs2e( SC, sclk_str);
% % 
% utc = cspice_et2utc( et, 'C', 6 );
% cspice_kclear;

%%
switch upper(crism_obs.info.obs_classType)
    case {'FFC'}
        line_idxes = 100:200;
        % line_idxes = 634:670;
    case {'FRT','HRL','FRS','HRS','ATO'}
        TRRRAdata = CRISMdata(basenameRA,'');
        [valid_lines,valid_samples] = crism_examine_valid_LinesColumns(TRRRAdata);
        valid_lines = find(valid_lines);
        line_idxes = valid_lines;
    otherwise
        error('Not implemented for %s',upper(crism_obs.info.obs_classType));
end

% x = 623; y = 432; x = 40; y = 229;
TRRIFdata.load_basenamesCDR();
DMdata = TRRIFdata.readCDR('DM');
DMdata.readimgi();
dm_mask_1nan = convertBoolTo1nan(DMdata.img==1);
[v_min,i_min] = min2d(DEdata.ddr.Elevation.img.*dm_mask_1nan(1,:,150));
x_min = i_min(2); y_min = i_min(1);
lat_min = DEdata.ddr.Latitude.img(y_min,x_min);
lon_min = DEdata.ddr.Longitude.img(y_min,x_min);

[v_max,i_max] = max2d(DEdata.ddr.Elevation.img.*dm_mask_1nan(1,:,150));
x_max = i_max(2); y_max = i_max(1);
lat_max = DEdata.ddr.Latitude.img(y_max,x_max);
lon_max = DEdata.ddr.Longitude.img(y_max,x_max);
% ina_areoid = DEdata.ddr.INA_at_areoid.img(y,x);
% ema_areoid = DEdata.ddr.EMA_at_areoid.img(y,x);

lat_mean = nanmean(DEdata.ddr.Latitude.img(line_idxes,:).*dm_mask_1nan(1,:,150),'all');
lon_mean = nanmean(DEdata.ddr.Longitude.img(line_idxes,:).*dm_mask_1nan(1,:,150),'all');
ina_areoid = nanmean(DEdata.ddr.INA_at_areoid.img(line_idxes,:).*dm_mask_1nan(1,:,150),'all');
ema_areoid = nanmean(DEdata.ddr.EMA_at_areoid.img(line_idxes,:).*dm_mask_1nan(1,:,150),'all');

switch upper(crism_obs.info.obs_classType)
    case {'FFC'}
        lat = lat_mean;
        lon = lon_mean;
    case {'FRT','HRL','FRS','HRS','ATO'}
        lat = lat_mean;
        lon = lon_mean;
    otherwise
        error('Not implemented for %s', ...
            upper(crism_obs.info.obs_classType));
end

lat = lat_mean;
lon = lon_mean;

[~,i2d] = min2d((DEdata.ddr.Latitude.img-lat).^2 + (DEdata.ddr.Longitude.img-lon).^2);

ema_areoid_ctr = DEdata.ddr.EMA_at_areoid.img(i2d(1),i2d(2));

if verbose
    TRRIFdata.set_rgb();
    X = DEdata.ddr.Longitude.img; Y = DEdata.ddr.Latitude.img;
    Z = DEdata.ddr.Elevation.img;
    figure;
    surf(X,Y,Z,TRRIFdata.RGB.CData_Scaled,'EdgeColor','none');
    % ax = gca;
end


%%
[pList] = design_presList(dt, lon, lat, ...
    'STEP',-30,'P_TOL',0.1,'TOL_M',1.05, 'scena', scena, 'MCD_VER', mcd_ver);
[params] = mcd_get_vertical_profile(dt,lon,lat,...
    'XZ',pList,'zkey',4,'verbose',false,'scena', scena, 'MCD_VER', mcd_ver);
% [pListMY30] = design_presList(dt,lon,lat,'STEP',-30,'P_TOL',0.1,'TOL_M',1.05,'scena',28);
% [paramsMY30] = mcd_get_vertical_profile(dt,lon,lat,...
%     'XZ',pListMY30,'zkey',4,'verbose',false,'scena',28);
alt_surf = [params.alt_above_surf];
alt_surf_mid = (alt_surf(2:end) + alt_surf(1:(end-1)))/2;
[params2] = mcd_get_vertical_profile(dt,lon,lat,...
    'XZ',alt_surf_mid,'zkey',3,'verbose', verbose, 'scena', scena, ...
    'MCD_VER', mcd_ver);

L_mid = diff(alt_surf);
T_mid = [params2.temp];
p_Pa_mid = [params2.pres];
co2_r_mid = [params2.vol_mix_ratio_co2];
co_r_mid = [params2.vol_mix_ratio_co];
h2ovapor_r_mid = [params2.vol_mix_ratio_h2ovapor];
alt_areoid_mid = [params2.alt_above_areoid];

R = 8.314462618;

n_mol = p_Pa_mid ./ (R.*T_mid) .* (L_mid*100.*10^(-4));

% 44.009756876395
% 28.0102330819987664
% 18.015266585696754


col_co2 = sum(44.009756876395*(n_mol.*co2_r_mid));
fprintf('CO2 column kg/m^2 = %f (database: %f).\n', ...
    col_co2, params2(1).col_co2);

col_h2o = sum(28.0102330819987664*(n_mol.*h2ovapor_r_mid));
fprintf('H2O vapor column kg/m^2 = %f.\n', col_h2o);

alt_areoid = [params.alt_above_areoid];

%% save
switch upper(crism_obs.info.obs_classType)
    case {'FFC'}
        fname = sprintf('%s_l%dt%d_MCDv%s_scena%d.mat', ...
            TRRIFdata.basename, line_idxes(1), line_idxes(end), ...
            mcd_ver, scena);
    case {'FRT','HRL','FRS','HRS','ATO'}
        fname = sprintf('%s_mean_MCDv%s_scena%d.mat', ...
            TRRIFdata.basename, mcd_ver, scena);
    otherwise
        error('Not implemented for %s', ...
            upper(crism_obs.info.obs_classType));
end

if save_file
    fpath = fullfile(save_dir, fname);
    save(fpath,'alt_surf', 'alt_surf_mid', 'alt_areoid', 'L_mid', ...
        'T_mid','p_Pa_mid', 'co2_r_mid', 'co_r_mid', 'h2ovapor_r_mid', ...
        'ina_areoid', 'ema_areoid', 'ema_areoid_ctr', 'scena', ...
        'mcd_ver', '-v7');
end

out = [];
out.alt_surf = alt_surf;
out.alt_surf_mid = alt_surf_mid;
out.alt_areoid = alt_areoid;
out.L_mid = L_mid;
out.T_mid = T_mid;
out.p_Pa_mid = p_Pa_mid;
out.co2_r_mid = co2_r_mid;
out.co_r_mid = co_r_mid;
out.h2ovapor_r_mid = h2ovapor_r_mid;
out.ina_areoid = ina_areoid;
out.ema_areoid = ema_areoid;
out.ema_areoid_ctr = ema_areoid_ctr;
out.mcd_ver = mcd_ver;
out.scena = scena;

[date_julian] = mcd_get_julian_date(dt,'VERBOSE', verbose, 'MCD_VER', mcd_ver);
out.date_julian = date_julian;

if save_file
    out.fpath = fpath;
end

end