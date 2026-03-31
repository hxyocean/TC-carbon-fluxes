%% this program is an example to correct the ccmp wind dataset with TC wind.
tic;
clear;

%% init
path_ccmp = './200007';
path_save = './';

climate_date_str = '20000721'; % the date of ccmp used to correction
climate_date = datenum(climate_date_str, 'yyyymmdd');

lon_ccmp = ncread([path_ccmp, '/CCMP_Wind_Analysis_', climate_date_str, '_V02.0_L3.0_RSS.nc'], 'longitude');
lat_ccmp = ncread([path_ccmp, '/CCMP_Wind_Analysis_', climate_date_str, '_V02.0_L3.0_RSS.nc'], 'latitude');
[grid_lat_ccmp, grid_lon_ccmp] = meshgrid(lat_ccmp, lon_ccmp);
grid_lat_ccmp = cast(grid_lat_ccmp, 'double');
grid_lon_ccmp = cast(grid_lon_ccmp, 'double');

%==========================================================================%
% read the TC best track data 

occur_TC_id       = zeros(4, 20);
occur_TC_idid     = zeros(4, 20);
occur_TC_count    = zeros(4, 1);

besttrack_path='../IBTrACS/IBTrACS.ALL.v04r00_latest.nc';
besttrack_info=ncinfo(besttrack_path);
besttrack_lat=ncread(besttrack_path,'lat');
besttrack_lon=ncread(besttrack_path,'lon');
besttrack_pre=ncread(besttrack_path,'wmo_pres');
besttrack_intens=ncread(besttrack_path,'wmo_wind');
besttrack_transspd=ncread(besttrack_path,'storm_speed');
besttrack_transdir=ncread(besttrack_path,'storm_dir');
besttrack_landfall=ncread(besttrack_path,'landfall');
besttrack_isotime=ncread(besttrack_path,'iso_time');
besttrack_dist2land=ncread(besttrack_path,'dist2land');
besttrack_lon(besttrack_lon<0)=besttrack_lon(besttrack_lon<0) + 360; % tansfer the range of longintude from [-180,180] to [0,360)
besttrack_intens = besttrack_intens * 1.852 / 3.6;                   % change the units from knot to m/s
besttrack_transspd = besttrack_transspd * 1.852 / 3.6;

besttrack_intens_back = besttrack_intens;

for i = 1:size(besttrack_intens, 2)    
    tem_num = find(isnan(besttrack_lat(:, i)), 1, "first") - 1;
    if datenum(besttrack_isotime(:, tem_num, i)') < climate_date ||...
            datenum(besttrack_isotime(:, 1, i)') >= climate_date + 1
        continue;
    end
    
    tem_isotime = datenum(besttrack_isotime(:, 1:tem_num, i)');
    for j = 1:tem_num
        if tem_isotime(j) < climate_date || tem_isotime(j) >= climate_date + 1 ||...
                mod((tem_isotime(j) - climate_date) * 4 + 1, 1) ~= 0
            continue;
        end
        occur_TC_count((tem_isotime(j) - climate_date) * 4 + 1) = ...
            occur_TC_count((tem_isotime(j) - climate_date) * 4 + 1) + 1;
        tem_count = occur_TC_count((tem_isotime(j) - climate_date) * 4 + 1);
        occur_TC_id((tem_isotime(j) - climate_date) * 4 + 1, tem_count) = i;
        occur_TC_idid((tem_isotime(j) - climate_date) * 4 + 1, tem_count) = j;
    end
end
%%  

lim_dist_out = 1200;
lim_dist_in  = 600;
ave_num = 11; % must be odd number
tem_store_u = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4 * ave_num);
tem_store_v = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4 * ave_num);

compare_wind = [];
%==========================================================================%
% init the store data
for i = 1:ave_num
    tem_current_time=datestr(climate_date + i - 1 - floor(ave_num / 2), 'yyyymmdd');
    try
        tem_u = ncread(fullfile(path_ccmp,...
            ['CCMP_Wind_Analysis_', tem_current_time, '_V02.0_L3.0_RSS.nc']), 'uwnd');
        tem_v = ncread(fullfile(path_ccmp,...
            ['CCMP_Wind_Analysis_', tem_current_time, '_V02.0_L3.0_RSS.nc']), 'vwnd');
    catch
            disp(['data missing: ', tem_current_time])
    end
        
    tem_store_u(:, :, (i-1)*4 + 1:i*4) = tem_u;
    tem_store_v(:, :, (i-1)*4 + 1:i*4) = tem_v;
end


%% calculate

u_filt = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4);
v_filt = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4);
u_tcwind = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4);
v_tcwind = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4);
for j = 1:4

    tem_mean_u = mean(tem_store_u(:, :, j:40+j), 3);
    tem_mean_v = mean(tem_store_v(:, :, j:40+j), 3);

    tem_mask = zeros(size(grid_lon_ccmp));
    tem_w    = zeros(size(grid_lon_ccmp));
    tem_TC_u = zeros(size(grid_lon_ccmp));
    tem_TC_v = zeros(size(grid_lon_ccmp));
    tem_maxwind = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), occur_TC_count(j));

    for k = 1:occur_TC_count(j)
        tem_lat = besttrack_lat(occur_TC_idid(j, k), occur_TC_id(j, k));
        tem_lon = besttrack_lon(occur_TC_idid(j, k), occur_TC_id(j, k));
        tem_vmax = besttrack_intens(occur_TC_idid(j, k), occur_TC_id(j, k));
        tem_maxwind_back = zeros(size(grid_lon_ccmp));

        
        % build TC wind
        if isnan(tem_vmax)
            % there is no record for TC maximum wind speed in IBTrACS dataset;
            % directly use the ccmp wind;
            u_wind = tem_store_u(:, :, j+20) - tem_mean_u;
            v_wind = tem_store_v(:, :, j+20) - tem_mean_v;
            tem_dist = fun_TC_dist(tem_lat, tem_lon, grid_lat_ccmp, grid_lon_ccmp); 
            tem_TC_u(tem_dist <= lim_dist_out) = ...
                tem_TC_u(tem_dist <= lim_dist_out) + u_wind(tem_dist <= lim_dist_out);
            tem_TC_v(tem_dist <= lim_dist_out) = ...
                tem_TC_v(tem_dist <= lim_dist_out) + v_wind(tem_dist <= lim_dist_out);
            tem_mask(tem_dist <= lim_dist_out) = 1;
            tem_w(tem_dist <= lim_dist_in) = 1;
            tem_maxwind_back(tem_dist <= lim_dist_in) = k;
        else
            [u_wind, v_wind, tem_dist] = ...
                fun_TC_wind(grid_lon_ccmp, grid_lat_ccmp, ...
                tem_lon, tem_lat, tem_vmax);
            tem_TC_u(tem_dist <= lim_dist_out) = ...
                tem_TC_u(tem_dist <= lim_dist_out) + u_wind(tem_dist <= lim_dist_out);
            tem_TC_v(tem_dist <= lim_dist_out) = ...
                tem_TC_v(tem_dist <= lim_dist_out) + v_wind(tem_dist <= lim_dist_out);
            tem_mask(tem_dist <= lim_dist_out) = 1;
            tem_w(tem_dist <= lim_dist_in) = 1;
            tem_maxwind_back(tem_dist <= lim_dist_in) = k;
        end
        tem_maxwind(:, :, k) = tem_maxwind_back;

    end
    
    tem_index_nan = (tem_w == 0 & tem_mask == 1);
    tem_index = ~ tem_index_nan;
    tem_lat = grid_lat_ccmp(tem_index);
    tem_lon = grid_lon_ccmp(tem_index);
    tem_data = tem_w(tem_index);
    
    F = scatteredInterpolant(tem_lat(:), tem_lon(:), tem_data(:), 'linear');

    tem_lat = grid_lat_ccmp(tem_index_nan);
    tem_lon = grid_lon_ccmp(tem_index_nan);
    tem_w(tem_index_nan) = F(tem_lat(:), tem_lon(:));
    
    % add TC wind


    u_filt(:, :, j) = tem_mean_u .* tem_w + tem_store_u(:, :, j+20) .* (1-tem_w);
    v_filt(:, :, j) = tem_mean_v .* tem_w + tem_store_v(:, :, j+20) .* (1-tem_w);
    
    
    tem_mean_u = tem_mean_u + tem_TC_u;
    tem_mean_v = tem_mean_v + tem_TC_v;
    
    u_tcwind(:, :, j) = tem_mean_u .* tem_w + tem_store_u(:, :, j+20) .* (1-tem_w);
    v_tcwind(:, :, j) = tem_mean_v .* tem_w + tem_store_v(:, :, j+20) .* (1-tem_w);


    %==================================================================%
    % record the information
    tem_wind_ccmp = sqrt(tem_store_u(:, :, j+20) .* tem_store_u(:, :, j+20) + ...
        tem_store_v(:, :, j+20) .* tem_store_v(:, :, j+20));
    tem_wind_tc   = sqrt(u_tcwind(:, :, j) .* u_tcwind(:, :, j) + ...
        v_tcwind(:, :, j) .* v_tcwind(:, :, j));   
    tem_wind_filt = sqrt(u_filt(:, :, j) .* u_filt(:, :, j) + ...
        v_filt(:, :, j) .* v_filt(:, :, j));
    for k = 1:occur_TC_count(j)
        tem_vmax = besttrack_intens(occur_TC_idid(j, k), occur_TC_id(j, k));
        tem_index_maxwind = (tem_maxwind(:, :, k) == k);
        compare_wind = [compare_wind; [tem_vmax, ...
            max(tem_wind_ccmp(tem_index_maxwind), [], 'all'), ...
            max(tem_wind_tc(tem_index_maxwind), [], 'all'), ...
            max(tem_wind_filt(tem_index_maxwind), [], 'all')]];
    end
end

%% save the wind data after correction
fun_save_wind(path_save, climate_date_str, ...
    sqrt(mean(u_tcwind.^2 + v_tcwind.^2, 3)), ...
    sqrt(mean(u_filt.^2 + v_filt.^2, 3)), ... 
    sqrt(mean(tem_store_u(:, :, 21: 24) .^ 2 + tem_store_v(:, :, 21: 24) .^ 2, 3)), ...
    lon_ccmp, lat_ccmp);

disp([climate_date_str, ' complete!'])
save([path_save, 'compare_wind.mat'], "compare_wind");

toc;