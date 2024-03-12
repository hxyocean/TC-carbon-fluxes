%% this program is used to calculated the fco2 induced by TC
tic;
clear;

addpath(genpath('../function'));
%% init
path_fco2_save = 'I:/ncs/data/data_rec/fco2_TC/';

path_pco2air = 'I:/ncs/data/data_rec/era_pressure/';
path_wind = 'I:/ncs/data/data_rec/ccmp_tcwind/';
path_pco2sea_rec = 'I:/ncs/data/data_rec/pco2sea_avhrr_hycom_smooth/';
path_sst = 'H:/AVHRR/AVHRR.v2.1/';
path_sss = 'I:/ncs/data/data_rec/hycom_sss_025/';
path_inform = 'I:/ncs/data/data_rec/tc_inform/all/';

%  grid
lon_sst = ncread('H:/AVHRR/AVHRR.v2.1/2000/oisst-avhrr-v02r01.20000101.nc', 'lon');
lat_sst = ncread('H:/AVHRR/AVHRR.v2.1/2000/oisst-avhrr-v02r01.20000101.nc', 'lat');
lon_sst = cast(lon_sst, 'double');
lat_sst = cast(lat_sst, 'double');
[grid_lat_sst, grid_lon_sst] = meshgrid(lat_sst, lon_sst);

% lon_pres = ncread('H:/era5 surface_pressure/1980/surface_pressure_19800101.nc', 'longitude');
% lat_pres = ncread('H:/era5 surface_pressure/1980/surface_pressure_19800101.nc', 'latitude');
% [grid_lat_pres, grid_lon_pres] = meshgrid(lat_pres, lon_pres);

lon_ccmp = ncread('H:/ccmp/ccmp/Y2000/CCMP_Wind_Analysis_20000101_V02.0_L3.0_RSS.nc', 'longitude');
lat_ccmp = ncread('H:/ccmp/ccmp/Y2000/CCMP_Wind_Analysis_20000101_V02.0_L3.0_RSS.nc', 'latitude');
lon_ccmp = cast(lon_ccmp, 'double');
lat_ccmp = cast(lat_ccmp, 'double');
[grid_lat_ccmp, grid_lon_ccmp] = meshgrid(lat_ccmp,lon_ccmp);

% basic grid
lon = lon_sst;
lat = lat_sst;
grid_lon = grid_lon_sst;
grid_lat = grid_lat_sst;

climate_timebegin = '1993-01-01';
climate_timeend   = '2020-12-31';

tem_num_year = str2double(climate_timeend(1:4)) - str2double(climate_timebegin(1:4)) + 1;

%================================================================================================%
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

% range_wind = nan(size(besttrack_intens, 2), 2);
% range_pre  = nan(size(besttrack_intens, 2), 2);
% for i = 1:size(besttrack_intens, 2)    
%     tem_num = find(isnan(besttrack_lat(:, i)), 1, "first") - 1;
%     if datenum(besttrack_isotime(:, tem_num, i)') < datenum(climate_timebegin) ||...
%             datenum(besttrack_isotime(:, 1, i)') >= datenum(climate_timeend) + 1
%         continue;
%     end
%     
%     tem_index_nan = (isnan(besttrack_intens(1:tem_num, i)));
%     tem_index_notnan = (~ tem_index_nan);
%     tem_index_nan = find(tem_index_nan);
%     tem_index_notnan = find(tem_index_notnan);
%     if length(tem_index_notnan) >= 2
%         tem_index = 1:tem_num;
%         besttrack_intens(tem_index_nan, i) = interp1(tem_index(tem_index_notnan), ...
%             besttrack_intens(tem_index_notnan, i)', tem_index(tem_index_nan), "linear");
%         besttrack_intens(besttrack_intens(:, i) < 0, i) = 0;
%     end
%     range_wind(i, 1) = min(besttrack_intens(:, i));
%     range_wind(i, 2) = max(besttrack_intens(:, i));
% 
%     %==========================================================================================%
%     tem_index_nan = (isnan(besttrack_pre(1:tem_num, i)));
%     tem_index_notnan = (~ tem_index_nan);
%     tem_index_nan = find(tem_index_nan);
%     tem_index_notnan = find(tem_index_notnan);
%     if length(tem_index_notnan) >= 2
%         tem_index = 1:tem_num;
%         besttrack_pre(tem_index_nan, i) = interp1(tem_index(tem_index_notnan), ...
%             besttrack_pre(tem_index_notnan, i)', tem_index(tem_index_nan), "linear");
%         besttrack_pre(besttrack_pre(:, i) < 0, i) = 0;
%     end
%     range_pre(i, 1) = min(besttrack_pre(:, i));
%     range_pre(i, 2) = max(besttrack_pre(:, i));
% end

climatology_sst = load('../../data/data_rec/climatology/climatology_sst.mat');
climatology_sst = climatology_sst.climate_sst;
climatology_sss = load('../../data/data_rec/climatology/climatology_sss.mat');
climatology_sss = climatology_sss.climate_sss;
climatology_pco2sea = load('../../data/data_rec/climatology/climatology_pco2sea.mat');
climatology_pco2sea = climatology_pco2sea.climate_pco2sea;
%==============================================================================================%
% loop year
loop_year = zeros(2023,1);
for i = 1:2023
    if (mod(i, 4) ==0 && mod(i, 100) ~= 0) || mod(i, 400) == 0
        loop_year(i) = 366;
    else
        loop_year(i) = 365;
    end
end
%==============================================================================================%
% area
m_area = zeros(size(grid_lon));
r_e = 6371 * 1000;
for i = 1:size(grid_lon, 1)
    for j = 1:size(grid_lon, 2)
        if grid_lat(i, j) < -40 || grid_lat(i, j) > 40
            continue;
        end
        m_area(i, j) = (r_e*cosd(grid_lat(i, j) - 0.125) + r_e*cosd(grid_lat(i, j) + 0.125)) * ...
            (0.25 / 180 * 3.1415926525) .* ...
            (r_e * 0.25 / 180 * 3.1415926525) / 2;
    end
end
%% main

limit_lat = 90;
gap_lat = 46;
limit_gap = 60;
before_time = 30;
rela_time = [-10, -4];
for i = 1:tem_num_year
    tem_current_year  = num2str(str2double(climate_timebegin(1:4)) + i - 1);
    tem_next_year     = num2str(str2double(climate_timebegin(1:4)) + i);
    tem_current_begin = datenum([tem_current_year, '-01-01']);
    tem_num_day_current = datenum([tem_next_year, '-01-01']) - datenum([tem_current_year, '-01-01']);
    tem_num_day = 365;

    tem_wind_tc      = zeros(size(grid_lon, 1), size(grid_lon, 2), tem_num_day + before_time);
    tem_wind_filt    = zeros(size(grid_lon, 1), size(grid_lon, 2), tem_num_day + before_time);
    tem_pco2air_tc   = zeros(size(grid_lon, 1), size(grid_lon, 2), tem_num_day + before_time);
    tem_pco2air_filt = zeros(size(grid_lon, 1), size(grid_lon, 2), tem_num_day + before_time);
    tem_pco2sea_tc      = zeros(size(grid_lon, 1), size(grid_lon, 2), tem_num_day + before_time);
    tem_pco2sea_cli  = zeros(size(grid_lon, 1), size(grid_lon, 2), tem_num_day + before_time);
    tem_SST_tc          = zeros(size(grid_lon, 1), size(grid_lon, 2), tem_num_day + before_time);
    tem_SSS_tc          = zeros(size(grid_lon, 1), size(grid_lon, 2), tem_num_day + before_time);
    tem_SST_cli      = zeros(size(grid_lon, 1), size(grid_lon, 2), tem_num_day + before_time);
    tem_SSS_cli      = zeros(size(grid_lon, 1), size(grid_lon, 2), tem_num_day + before_time);
    

    disp(['Find TC day:  ', tem_current_year]);
    tem_tc_inform = load([path_inform, tem_current_year, '.mat']);
    tem_tc_inform = tem_tc_inform.tem_tc_inform;


    disp(['Initiating data:  ', tem_current_year]);
    for j = (- before_time + 1) : tem_num_day
        tem_current_time = datestr(tem_current_begin + j - 1, 'yyyymmdd');
        if loop_year(str2double(tem_current_year)) == 366
            if j >= 60
                tem_current_time = datestr(tem_current_begin + j, 'yyyymmdd');
            else
                tem_current_time = datestr(tem_current_begin + j - 1, 'yyyymmdd');
            end
        end
        tem_wind_tc(:, gap_lat+1 : end-gap_lat, j + before_time)   = ncread([path_wind, ...
                tem_current_time(1:4), '/', 'CCMP_', tem_current_time, ...
                '_V2.0_TCWIND_REC.nc'], 'wind_tc_s');
        tem_wind_filt(:, gap_lat+1:end-gap_lat, j + before_time) = ncread([path_wind, ...
                tem_current_time(1:4), '/', 'CCMP_',tem_current_time, ...
                '_V2.0_TCWIND_REC.nc'], 'wind_filt');
        tem_data = load([path_pco2sea_rec, tem_current_time(1:4), '/', ...
            tem_current_time, '_pco2sea.mat']);
        tem_pco2sea_tc(:, :, j + before_time) = tem_data.pco2_sea;
        tem_pco2air_tc(:, :, j + before_time) = ncread([path_pco2air, ...
                tem_current_time(1:4), '/', 'surface_pressure_',tem_current_time, ...
                '_TCPRES_REC.nc'], 'pco2_tc');
        tem_pco2air_filt(:, :, j + before_time) = ncread([path_pco2air, ...
                tem_current_time(1:4), '/', 'surface_pressure_',tem_current_time, ...
                '_TCPRES_REC.nc'], 'pco2_filt');
        tem_SST_tc(:, :, j + before_time) = ncread([path_sst, ...
                tem_current_time(1:4), '/', 'oisst-avhrr-v02r01.',tem_current_time, ...
                '.nc'], 'sst');
        tem_data = load([path_sss, tem_current_time(1:4), '/', ...
            tem_current_time, '_sss.mat']);
        tem_SSS_tc(:, :, j + before_time) = tem_data.sss;
    end
    
    tem_SST_cli = cat(3, climatology_sst(:, :, end-before_time+1:end), ...
        climatology_sst);
    tem_SSS_cli = cat(3, climatology_sss(:, :, end-before_time+1:end), ...
        climatology_sss);
    tem_pco2sea_cli = cat(3, climatology_pco2sea(:, :, end-before_time+1:end), ...
        climatology_pco2sea);
    %=====================================================================%
    disp(['calculting:  ', tem_current_year]);
    fco2_wind = zeros(size(grid_lon, 1), size(grid_lon, 2));
    fco2_cool = zeros(size(grid_lon, 1), size(grid_lon, 2));
    fco2_all = zeros(size(grid_lon, 1), size(grid_lon, 2));
    for j = 1:size(grid_lon, 1)
        for k = 1:size(grid_lon, 2)
            tem_TC_time = find(tem_tc_inform(j, k, :));
            if isempty(tem_TC_time) || abs(lat_sst(k)) > limit_lat  
                continue
            end
            tem_TC_time_gap = tem_TC_time(2:end) - tem_TC_time(1:end-1);
            tem_TC_time_index = find(tem_TC_time_gap > limit_gap);
            tem_TC_time_begin = [tem_TC_time(1); tem_TC_time(tem_TC_time_index + 1)];
            tem_TC_time_end   = [tem_TC_time(tem_TC_time_index) + limit_gap; ...
                min(tem_TC_time(end) + limit_gap, tem_num_day)];
            
            for m = 1:length(tem_TC_time_begin)
                
                tem_b = tem_TC_time_begin(m) + before_time - before_time;
                tem_e = tem_TC_time_end(m) + before_time;

                tem_wind_tc_poi = tem_wind_tc(j, k, tem_b:tem_e);
                tem_wind_poi = tem_wind_filt(j, k, tem_b:tem_e);
                tem_pco2sea_tc_poi  = tem_pco2sea_tc(j, k, tem_b:tem_e);
                tem_pco2air_poi  = tem_pco2air_filt(j, k, tem_b:tem_e);
                tem_pco2air_tc_poi = tem_pco2air_tc(j, k, tem_b:tem_e);
                tem_SST_tc_poi = tem_SST_tc(j, k, tem_b:tem_e);
                tem_SSS_tc_poi = tem_SSS_tc(j, k, tem_b:tem_e);
                
                tem_pco2sea_poi = tem_pco2sea_cli(j, k, tem_b:tem_e);
                tem_pco2sea_poi = tem_pco2sea_poi + ...
                    mean(tem_pco2sea_tc_poi(before_time + rela_time(1) + 1:before_time + rela_time(2) + 1) - ...
                    tem_pco2sea_poi(before_time + rela_time(1) + 1:before_time + rela_time(2) + 1));
                tem_SST_poi = tem_SST_cli(j, k, tem_b:tem_e);
                tem_SST_poi = tem_SST_poi + ...
                    mean(tem_SST_tc_poi(before_time + rela_time(1) + 1:before_time + rela_time(2) + 1) - ...
                    tem_SST_poi(before_time + rela_time(1) + 1:before_time + rela_time(2) + 1));
                tem_SSS_poi = tem_SSS_cli(j, k, tem_b:tem_e);
                tem_SSS_poi = tem_SSS_poi + ...
                    mean(tem_SSS_tc_poi(before_time + rela_time(1) + 1:before_time + rela_time(2) + 1) - ...
                    tem_SSS_poi(before_time + rela_time(1) + 1:before_time + rela_time(2) + 1));
                
                % fco2
                tem_dfco2 = fun_fco2(tem_wind_tc_poi .* tem_wind_tc_poi, tem_pco2sea_tc_poi - tem_pco2air_tc_poi, ...
                    tem_SST_tc_poi, tem_SSS_tc_poi, "instant_14") - ...
                    fun_fco2(tem_wind_poi .* tem_wind_poi, tem_pco2sea_poi - tem_pco2air_poi, ...
                    tem_SST_poi, tem_SSS_poi, "instant_14");
                fco2_all(j, k) = fco2_all(j, k) + sum(squeeze(tem_dfco2(before_time + rela_time(2) + 2:end))) * 24;


                tem_dfco2 = fun_fco2(tem_wind_tc_poi .* tem_wind_tc_poi, tem_pco2sea_poi - tem_pco2air_tc_poi, ...
                    tem_SST_poi, tem_SSS_poi, "instant_14") - ...
                    fun_fco2(tem_wind_poi .* tem_wind_poi, tem_pco2sea_poi - tem_pco2air_poi, ...
                    tem_SST_poi, tem_SSS_poi, "instant_14");
                fco2_wind(j, k) = fco2_wind(j, k) + sum(squeeze(tem_dfco2(before_time + rela_time(2) + 2:end))) * 24;

                tem_dfco2 = fun_fco2(tem_wind_tc_poi .* tem_wind_tc_poi, tem_pco2sea_tc_poi - tem_pco2air_tc_poi, ...
                    tem_SST_tc_poi, tem_SSS_tc_poi, "instant_14") - ...
                    fun_fco2(tem_wind_tc_poi .* tem_wind_tc_poi, tem_pco2sea_poi - tem_pco2air_tc_poi, ...
                    tem_SST_poi, tem_SSS_poi, "instant_14");
                fco2_cool(j, k) = fco2_cool(j, k) + sum(squeeze(tem_dfco2(before_time + rela_time(2) + 2:end))) * 24;
            end

        end
    end
    
    fco2 = fun_fco2(tem_wind_filt .* tem_wind_filt, tem_pco2sea_tc - tem_pco2air_filt, ...
        tem_SST_tc, tem_SSS_tc, "instant_14") * 24;
    fco2 = fco2(:,:, before_time + 1 : end);
    fco2_season = cat(2, sum(fco2(:, 1:360, [1:120, 335:365]), 3), ...
        sum(fco2(:, 361:720, [152:304]), 3));
    fco2 = sum(fco2, 3);
    
%     parsave([path_fco2_save, tem_current_year, '.mat'], ...
%         fco2_cool, fco2_wind, fco2_all, fco2, fco2_season);
    sprintf('%d: %f %f %f %f %f', str2double(tem_current_year), ...
        nansum(m_area.*fco2_all, 'all') * 12 / 1e15, ...
        nansum(m_area.*fco2_wind, 'all') * 12 / 1e15, ...
        nansum(m_area.*fco2_cool, 'all') * 12 / 1e15, ...
        nansum(m_area.*fco2, 'all') * 12 / 1e15, ...
        nansum(m_area.*fco2_season, 'all') * 12 / 1e15)
end

toc;
