%% this program is used to store the point within the 200 km of TC center
tic;
clear;
addpath(genpath('../function'));
%% init
path_save = 'I:/ncs/data/data_rec/tc_inform/all/25/';

lon_sst = ncread('H:/AVHRR/AVHRR.v2.1/2000/oisst-avhrr-v02r01.20000101.nc', 'lon');
lat_sst = ncread('H:/AVHRR/AVHRR.v2.1/2000/oisst-avhrr-v02r01.20000101.nc', 'lat');
lon_sst = cast(lon_sst, 'double');
lat_sst = cast(lat_sst, 'double');
[grid_lat_sst, grid_lon_sst] = meshgrid(lat_sst, lon_sst);

% best track data
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


grid_lon = grid_lon_sst;
grid_lat = grid_lat_sst;
climate_timebegin = '1992-01-01';
climate_timeend   = '2021-12-31';
tem_num_year = str2double(climate_timeend(1:4)) - str2double(climate_timebegin(1:4)) + 1;
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
%% calculate
% 
lim_dist = 25;
delete(gcp('nocreate'));         % stop the process before start a new run.
numCore = feature('numcores');   % get the maxmium core num of PC.
parpool(numCore-13);                % start parpool.

parfor i = 1:tem_num_year

    tem_current_year  = num2str(str2double(climate_timebegin(1:4)) + i - 1);
    tem_next_year     = num2str(str2double(climate_timebegin(1:4)) + i);
    tem_time_begin    = datenum([tem_current_year, '-01-01']);
    tem_time_end      = datenum([tem_next_year, '-01-01']);
    tem_num_day = tem_time_end - tem_time_begin;
    tem_num_day = 365;
    tem_tc_inform          = zeros(size(grid_lon, 1), size(grid_lon, 2), tem_num_day, 'logical');
    tem_tc_inform_num      = zeros(size(grid_lon, 1), size(grid_lon, 2), tem_num_day, 'int16');
    tem_tc_inform_index    = zeros(size(grid_lon, 1), size(grid_lon, 2), tem_num_day, 3, 'int16');
    
    for j = 1:size(besttrack_intens, 2)
        tem_num = find(isnan(besttrack_lat(:, j)), 1, "first") - 1;
        
        if datenum(besttrack_isotime(:, tem_num, j)') < tem_time_begin ||...
                datenum(besttrack_isotime(:, 1, j)') >= tem_time_end
            continue;
        end
        tem_isotime = datenum(besttrack_isotime(:, 1:tem_num, j)');
        for k = 1:tem_num
            if tem_isotime(k) < tem_time_begin || tem_isotime(k) >= tem_time_end
                continue;
            end
            tem_index = floor(tem_isotime(k) - tem_time_begin + 1);
            if loop_year(str2double(tem_current_year)) == 366
                if tem_index == 60
                    continue;
                elseif tem_index > 60
                    tem_index = tem_index - 1;
                end
            end
            
            tem_lat = besttrack_lat(k, j);
            tem_lon = besttrack_lon(k, j);
            tem_dist = fun_TC_dist(tem_lat, tem_lon, grid_lat_sst, grid_lon_sst);
%             tem_dist = 6371 .* acos(sind(tem_lat).*sind(grid_lat_sst)+...
%                     cosd(tem_lat).*cosd(grid_lat_sst).*cosd(grid_lon_sst-tem_lon));

            tem_data = tem_tc_inform(:, :, tem_index);
            tem_data(tem_dist <= lim_dist) = true;
            tem_tc_inform(:, :, tem_index) = tem_data;

            [tem_row, tem_col] = find(tem_dist <= lim_dist);
            for m = 1:length(tem_row)
                tem_num_tc = tem_tc_inform_num(tem_row(m), tem_col(m), tem_index);
                if tem_num_tc == 0
                    tem_num_tc = tem_num_tc + 1;
                    tem_tc_inform_num(tem_row(m), tem_col(m), tem_index) = tem_num_tc;
                    tem_tc_inform_index(tem_row(m), tem_col(m), tem_index, tem_num_tc) = j;
                else
                    if tem_tc_inform_index(tem_row(m), tem_col(m), tem_index, tem_num_tc) == j
                        continue;
                    end
                    tem_num_tc = tem_num_tc + 1;
                    tem_tc_inform_num(tem_row(m), tem_col(m), tem_index) = tem_num_tc;
                    tem_tc_inform_index(tem_row(m), tem_col(m), tem_index, tem_num_tc) = j;
                end
            end
        end
        
    end
    parsave_large([path_save, tem_current_year, '.mat'], tem_tc_inform, tem_tc_inform_num, tem_tc_inform_index);
    disp([tem_current_year, '  completed!']);
end

% tem_data = zeros(size(grid_lon_sst));
% for i = 1:tem_num_day
%     tem_data(tem_tc_inform(:, :, i)) = tem_data(tem_tc_inform(:, :, i)) + 1; 
% end
% contourf(grid_lon, grid_lat, tem_data, 0:20);
% colorbar();

toc;