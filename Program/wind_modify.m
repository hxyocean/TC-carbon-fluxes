%% this program is to add the TC wind to ccmp
tic;
clear;

%% init
path_save = 'H:/ccmp/TC_wind/';
path_save = 'I:/ncs/data/data_rec/ccmp_tcwind/';
% path_save = 'I:/ncs/data/data_rec/case/';

path_wind = 'H:/ccmp/';
lon_ccmp = ncread('H:/ccmp/ccmp/Y2000/CCMP_Wind_Analysis_20000101_V02.0_L3.0_RSS.nc', 'longitude');
lat_ccmp = ncread('H:/ccmp/ccmp/Y2000/CCMP_Wind_Analysis_20000101_V02.0_L3.0_RSS.nc', 'latitude');
[grid_lat_ccmp, grid_lon_ccmp] = meshgrid(lat_ccmp,lon_ccmp);
grid_lat_ccmp = cast(grid_lat_ccmp, 'double');
grid_lon_ccmp = cast(grid_lon_ccmp, 'double');


climate_timebegin=datenum('1995-08-12');
climate_timeend=datenum('1995-08-15');

%==========================================================================%
% read the TC best track data 
besttrack_path='../../data/data_statistics/IBTrACS/IBTrACS.ALL.v04r00_latest.nc';
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
besttrack_rmw = ncread(besttrack_path,'usa_rmw') * 1.852;


besttrack_lon(besttrack_lon<0)=besttrack_lon(besttrack_lon<0) + 360; % tansfer the range of longintude from [-180,180] to [0,360)
besttrack_intens=besttrack_intens * 1.852 / 3.6;
besttrack_transspd = besttrack_transspd * 1.852 / 3.6;

besttrack_intens_back = besttrack_intens;

data_num_days = climate_timeend - climate_timebegin + 1;
occur_TC_id       = zeros(data_num_days * 4, 20);
occur_TC_idid     = zeros(data_num_days * 4, 20);
occur_TC_count    = zeros(data_num_days * 4, 1);

range_wind = nan(size(besttrack_intens, 2), 2);
range_rmw  = nan(size(besttrack_intens, 2), 2);
for i = 1:size(besttrack_intens, 2)    
    tem_num = find(isnan(besttrack_lat(:, i)), 1, "first") - 1;
    if datenum(besttrack_isotime(:, tem_num, i)') < climate_timebegin ||...
            datenum(besttrack_isotime(:, 1, i)') >= climate_timeend + 1
        continue;
    end
    %==========================================================================================%
    tem_index_nan = (isnan(besttrack_intens(1:tem_num, i)));
    tem_index_notnan = (~ tem_index_nan);
    tem_index_nan = find(tem_index_nan);
    tem_index_notnan = find(tem_index_notnan);
    if length(tem_index_notnan) >= 2
        tem_index = 1:tem_num;
        besttrack_intens(tem_index_nan, i) = interp1(tem_index(tem_index_notnan), ...
            besttrack_intens(tem_index_notnan, i)', tem_index(tem_index_nan), "linear");
        besttrack_intens(besttrack_intens(:, i) < 0, i) = 0;
    end
    range_wind(i, 1) = min(besttrack_intens(:, i));
    range_wind(i, 2) = max(besttrack_intens(:, i));
    %==========================================================================================%
    tem_index_nan = (isnan(besttrack_rmw(1:tem_num, i)));
    tem_index_notnan = (~ tem_index_nan);
    tem_index_nan = find(tem_index_nan);
    tem_index_notnan = find(tem_index_notnan);
    if length(tem_index_notnan) >= 2
        tem_index = 1:tem_num;
        besttrack_rmw(tem_index_nan, i) = interp1(tem_index(tem_index_notnan), ...
            besttrack_rmw(tem_index_notnan, i)', tem_index(tem_index_nan), "linear");
        besttrack_rmw(besttrack_rmw(:, i) < 0, i) = 0;
    end
    range_rmw(i, 1) = min(besttrack_rmw(:, i));
    range_rmw(i, 2) = max(besttrack_rmw(:, i));
    %==========================================================================================%
    tem_isotime = datenum(besttrack_isotime(:, 1:tem_num, i)');
    for j = 1:tem_num
        if tem_isotime(j) < climate_timebegin || tem_isotime(j) >= climate_timeend+1 ||...
                mod((tem_isotime(j) - climate_timebegin) * 4 + 1, 1) ~= 0
            continue;
        end
        occur_TC_count((tem_isotime(j) - climate_timebegin) * 4 + 1) = ...
            occur_TC_count((tem_isotime(j) - climate_timebegin) * 4 + 1) + 1;
        tem_count = occur_TC_count((tem_isotime(j) - climate_timebegin) * 4 + 1);
        occur_TC_id((tem_isotime(j) - climate_timebegin) * 4 + 1, tem_count) = i;
        occur_TC_idid((tem_isotime(j) - climate_timebegin) * 4 + 1, tem_count) = j;
    end
end
%% 

lim_dist_out = 1200;
lim_dist_in  = 600;
lim_dist_maxwind = 400;
ave_num = 11; % must be odd number
tem_store_u = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4 * ave_num);
tem_store_v = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4 * ave_num);

compare_wind = [];
%==========================================================================%
% init the store data
for i = 1:ave_num
    tem_current_time=datestr(climate_timebegin + i - 1 - floor(ave_num / 2),'yyyymmdd');
    if str2double(tem_current_time(1:4)) < 2019
        tem_u = ncread(fullfile(path_wind,['ccmp/Y',tem_current_time(1:4)],...
            ['CCMP_Wind_Analysis_', tem_current_time, '_V02.0_L3.0_RSS.nc']), 'uwnd');
        tem_v = ncread(fullfile(path_wind,['ccmp/Y',tem_current_time(1:4)],...
            ['CCMP_Wind_Analysis_', tem_current_time, '_V02.0_L3.0_RSS.nc']), 'vwnd');
    else
        try
            tem_u = ncread(fullfile(path_wind,['ccmp_v2.1_NRT/Y',tem_current_time(1:4)],...
                ['CCMP_RT_Wind_Analysis_', tem_current_time, '_V02.1_L3.0_RSS.nc']), 'uwnd');
            tem_v = ncread(fullfile(path_wind,['ccmp_v2.1_NRT/Y',tem_current_time(1:4)],...
                ['CCMP_RT_Wind_Analysis_', tem_current_time, '_V02.1_L3.0_RSS.nc']), 'vwnd');
        catch
            if string(tem_current_time) == "20201022" || string(tem_current_time) == "20201023"
                tem_u_b = ncread(fullfile(path_wind,['ccmp_v2.1_NRT/Y',tem_current_time(1:4)],...
                    'CCMP_RT_Wind_Analysis_20201021_V02.1_L3.0_RSS.nc'), 'uwnd');
                tem_v_b = ncread(fullfile(path_wind,['ccmp_v2.1_NRT/Y',tem_current_time(1:4)],...
                    'CCMP_RT_Wind_Analysis_20201021_V02.1_L3.0_RSS.nc'), 'vwnd');
                tem_u_a = ncread(fullfile(path_wind,['ccmp_v2.1_NRT/Y',tem_current_time(1:4)],...
                    'CCMP_RT_Wind_Analysis_20201024_V02.1_L3.0_RSS.nc'), 'uwnd');
                tem_v_a = ncread(fullfile(path_wind,['ccmp_v2.1_NRT/Y',tem_current_time(1:4)],...
                    'CCMP_RT_Wind_Analysis_20201024_V02.1_L3.0_RSS.nc'), 'vwnd');
                tem_u = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4, "single");
                tem_v = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4, "single");
                for j = 1:4
                    tem_time =  datenum(tem_current_time,'yyyymmdd') - datenum('2020-10-22') + j * 0.25;
                    tem_u(:, :, j) = tem_u_b(:, :, end) + ...
                        (tem_u_a(:, :, 1) - tem_u_b(:, :, end)) ./ 2.25 .* tem_time;
                    tem_v(:, :, j) = tem_v_b(:, :, end) + ...
                        (tem_v_a(:, :, 1) - tem_v_b(:, :, end)) ./ 2.25 .* tem_time;
                end
            else
                tem_u = ncread(fullfile(path_wind,['ccmp_v2.0_NRT/Y',tem_current_time(1:4)],...
                    ['CCMP_RT_Wind_Analysis_', tem_current_time, '_V02.0_L3.0_RSS.nc']), 'uwnd');
                tem_v = ncread(fullfile(path_wind,['ccmp_v2.0_NRT/Y',tem_current_time(1:4)],...
                    ['CCMP_RT_Wind_Analysis_', tem_current_time, '_V02.0_L3.0_RSS.nc']), 'vwnd');
            end
        end
        
    end
    tem_store_u(:, :, (i-1)*4 + 1:i*4) = tem_u;
    tem_store_v(:, :, (i-1)*4 + 1:i*4) = tem_v;
end

%==========================================================================%
% calculate
for i = 1:(climate_timeend - climate_timebegin + 1)
    tem_current_time=datestr(climate_timebegin+i-1,'yyyymmdd');
    u_filt = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4);
    v_filt = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4);
    u_tcwind = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4);
    v_tcwind = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4);
    u_tcwind_s = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4);
    v_tcwind_s = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4);
    for j = 1:4

        tem_mean_u = mean(tem_store_u(:, :, j:40+j), 3);
        tem_mean_v = mean(tem_store_v(:, :, j:40+j), 3);

        tem_index = (i - 1) * 4 + j;
        tem_mask = zeros(size(grid_lon_ccmp));
        tem_w    = zeros(size(grid_lon_ccmp));
        tem_TC_u = zeros(size(grid_lon_ccmp));
        tem_TC_v = zeros(size(grid_lon_ccmp));
        tem_maxwind = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), occur_TC_count(tem_index));
        for k = 1:occur_TC_count(tem_index)
            tem_lat = besttrack_lat(occur_TC_idid(tem_index, k), occur_TC_id(tem_index, k));
            tem_lon = besttrack_lon(occur_TC_idid(tem_index, k), occur_TC_id(tem_index, k));
            tem_vmax = besttrack_intens(occur_TC_idid(tem_index, k), occur_TC_id(tem_index, k));
            tem_rmw = besttrack_rmw(occur_TC_idid(tem_index, k), occur_TC_id(tem_index, k));
            tem_maxwind_back = zeros(size(grid_lon_ccmp));
            if isnan(tem_vmax)
                u_wind = tem_store_u(:, :, j+20) - tem_mean_u;
                v_wind = tem_store_v(:, :, j+20) - tem_mean_v;
                tem_dist = fun_TC_dist(tem_lat, tem_lon, grid_lat_ccmp, grid_lon_ccmp); 
                tem_TC_u(tem_dist <= lim_dist_out) = ...
                    tem_TC_u(tem_dist <= lim_dist_out) + u_wind(tem_dist <= lim_dist_out);
                tem_TC_v(tem_dist <= lim_dist_out) = ...
                    tem_TC_v(tem_dist <= lim_dist_out) + v_wind(tem_dist <= lim_dist_out);
                tem_mask(tem_dist <= lim_dist_out) = 1;
                tem_w(tem_dist <= lim_dist_in) = 1;
                tem_maxwind_back(tem_dist <= lim_dist_maxwind) = k;
            else
                [u_wind, v_wind, tem_dist] = ...
                    fun_TC_wind(grid_lon_ccmp, grid_lat_ccmp, ...
                    tem_lon, tem_lat, tem_vmax, tem_rmw);
                tem_TC_u(tem_dist <= lim_dist_out) = ...
                    tem_TC_u(tem_dist <= lim_dist_out) + u_wind(tem_dist <= lim_dist_out);
                tem_TC_v(tem_dist <= lim_dist_out) = ...
                    tem_TC_v(tem_dist <= lim_dist_out) + v_wind(tem_dist <= lim_dist_out);
                tem_mask(tem_dist <= lim_dist_out) = 1;
                tem_w(tem_dist <= lim_dist_in) = 1;
                tem_maxwind_back(tem_dist <= lim_dist_maxwind) = k;
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

        tem_w(tem_mean_u.*tem_mean_u + tem_mean_v.*tem_mean_v < ...
            tem_store_u(:, :, j+20).*tem_store_u(:, :, j+20)+ ...
            tem_store_v(:, :, j+20).*tem_store_v(:, :, j+20)) = 0;
        u_tcwind_s(:, :, j) = tem_mean_u .* tem_w + tem_store_u(:, :, j+20) .* (1-tem_w);
        v_tcwind_s(:, :, j) = tem_mean_v .* tem_w + tem_store_v(:, :, j+20) .* (1-tem_w);

        %==================================================================%
        % record the information
        tem_wind_ccmp = sqrt(tem_store_u(:, :, j+20) .* tem_store_u(:, :, j+20) + ...
            tem_store_v(:, :, j+20) .* tem_store_v(:, :, j+20));
        tem_wind_tc_s = sqrt(u_tcwind_s(:, :, j) .* u_tcwind_s(:, :, j) + ...
            v_tcwind_s(:, :, j) .* v_tcwind_s(:, :, j));
        tem_wind_tc   = sqrt(u_tcwind(:, :, j) .* u_tcwind(:, :, j) + ...
            v_tcwind(:, :, j) .* v_tcwind(:, :, j));   
        tem_wind_filt = sqrt(u_filt(:, :, j) .* u_filt(:, :, j) + ...
            v_filt(:, :, j) .* v_filt(:, :, j));
        tem_index = (i - 1) * 4 + j;
        for k = 1:occur_TC_count(tem_index)
            tem_vmax = besttrack_intens(occur_TC_idid(tem_index, k), occur_TC_id(tem_index, k));
            tem_vmax_back = besttrack_intens_back(occur_TC_idid(tem_index, k), occur_TC_id(tem_index, k));
            tem_index_maxwind = (tem_maxwind(:, :, k) == k);
            compare_wind = [compare_wind; [tem_vmax_back, tem_vmax, ...
                max(tem_wind_ccmp(tem_index_maxwind), [], 'all'), ...
                max(tem_wind_tc(tem_index_maxwind), [], 'all'), ...
                max(tem_wind_tc_s(tem_index_maxwind), [], 'all'), ...
                max(tem_wind_filt(tem_index_maxwind), [], 'all')]];
        end
    end
%     fun_save(path_save, tem_current_time, u_tcwind_s, v_tcwind_s, u_filt, v_filt,...
%         sqrt(mean(u_tcwind_s.*u_tcwind_s + v_tcwind_s.*v_tcwind_s, 3)), ...
%         sqrt(mean(u_filt.*u_filt + v_filt.*v_filt, 3)), ... 
%         lon_ccmp, lat_ccmp);

%     fun_save_wind(path_save, tem_current_time, ...
%         sqrt(mean(u_tcwind.*u_tcwind + v_tcwind.*v_tcwind, 3)), ...
%         sqrt(mean(u_filt.*u_filt + v_filt.*v_filt, 3)), ... 
%         sqrt(mean(tem_store_u(:, :, 21: 24).*tem_store_u(:, :, 21: 24) + ...
%             tem_store_v(:, :, 21: 24).* tem_store_v(:, :, 21: 24), 3)), ...
%         sqrt(mean(u_tcwind_s.*u_tcwind_s + v_tcwind_s.*v_tcwind_s, 3)), ...
%         lon_ccmp, lat_ccmp);
    disp([tem_current_time, ' complete!'])
    %=======================================================================================%
    % update the data
    tem_current_time = datestr(climate_timebegin + i - 1 + floor(ave_num / 2) + 1,'yyyymmdd');
    if str2double(tem_current_time(1:4)) < 2019
        tem_u = ncread(fullfile(path_wind,['ccmp/Y',tem_current_time(1:4)],...
            ['CCMP_Wind_Analysis_', tem_current_time, '_V02.0_L3.0_RSS.nc']), 'uwnd');
        tem_v = ncread(fullfile(path_wind,['ccmp/Y',tem_current_time(1:4)],...
            ['CCMP_Wind_Analysis_', tem_current_time, '_V02.0_L3.0_RSS.nc']), 'vwnd');
    else
        try
            tem_u = ncread(fullfile(path_wind,['ccmp_v2.1_NRT/Y',tem_current_time(1:4)],...
                ['CCMP_RT_Wind_Analysis_', tem_current_time, '_V02.1_L3.0_RSS.nc']), 'uwnd');
            tem_v = ncread(fullfile(path_wind,['ccmp_v2.1_NRT/Y',tem_current_time(1:4)],...
                ['CCMP_RT_Wind_Analysis_', tem_current_time, '_V02.1_L3.0_RSS.nc']), 'vwnd');
        catch
            if string(tem_current_time) == "20201022" || string(tem_current_time) == "20201023"
                tem_u_b = ncread(fullfile(path_wind,['ccmp_v2.1_NRT/Y',tem_current_time(1:4)],...
                    'CCMP_RT_Wind_Analysis_20201021_V02.1_L3.0_RSS.nc'), 'uwnd');
                tem_v_b = ncread(fullfile(path_wind,['ccmp_v2.1_NRT/Y',tem_current_time(1:4)],...
                    'CCMP_RT_Wind_Analysis_20201021_V02.1_L3.0_RSS.nc'), 'vwnd');
                tem_u_a = ncread(fullfile(path_wind,['ccmp_v2.1_NRT/Y',tem_current_time(1:4)],...
                    'CCMP_RT_Wind_Analysis_20201024_V02.1_L3.0_RSS.nc'), 'uwnd');
                tem_v_a = ncread(fullfile(path_wind,['ccmp_v2.1_NRT/Y',tem_current_time(1:4)],...
                    'CCMP_RT_Wind_Analysis_20201024_V02.1_L3.0_RSS.nc'), 'vwnd');
                tem_u = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4, "single");
                tem_v = zeros(size(grid_lon_ccmp, 1), size(grid_lon_ccmp, 2), 4, "single");
                for j = 1:4
                    tem_time =  datenum(tem_current_time, 'yyyymmdd') - datenum('2020-10-22') + j * 0.25;
                    tem_u(:, :, j) = tem_u_b(:, :, end) + ...
                        (tem_u_a(:, :, 1) - tem_u_b(:, :, end)) ./ 2.25 .* tem_time;
                    tem_v(:, :, j) = tem_v_b(:, :, end) + ...
                        (tem_v_a(:, :, 1) - tem_v_b(:, :, end)) ./ 2.25 .* tem_time;
                end
            else
                tem_u = ncread(fullfile(path_wind,['ccmp_v2.0_NRT/Y',tem_current_time(1:4)],...
                    ['CCMP_RT_Wind_Analysis_', tem_current_time, '_V02.0_L3.0_RSS.nc']), 'uwnd');
                tem_v = ncread(fullfile(path_wind,['ccmp_v2.0_NRT/Y',tem_current_time(1:4)],...
                    ['CCMP_RT_Wind_Analysis_', tem_current_time, '_V02.0_L3.0_RSS.nc']), 'vwnd');
            end
        end
        
    end
    tem_store_u(:, :, 1:end-4) = tem_store_u(:, :, 5:end);
    tem_store_v(:, :, 1:end-4) = tem_store_v(:, :, 5:end);

    tem_store_u(:, :, end-3:end) = tem_u;
    tem_store_v(:, :, end-3:end) = tem_v;
    
end

save([path_save, 'compare_wind.mat'], "compare_wind");
toc;