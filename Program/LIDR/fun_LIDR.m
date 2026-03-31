%% this program is used to calculate the surface DIC value based on LIDR method
function [DIC] = fun_LIDR(sst, sss, lon, lat, messure_date)
% sst: sea surface temperature (Celsius)
% sss: sea surface salinity (psu)
% lon: longitude (degree)
% lat: latitude  (degree)
% messure_data: time of messurement ('yyyy-mm-dd')

% the input data must be single point or mesh grid

load('./interpolate_cor.mat', 'F0', 'F_SST', 'F_SSS');
load('./corr_modify.mat', "corr_b", "corr_mean_time");

%%
gap = 5.0;
lon = mod(lon, 360);

grid_F0 = F0(lon, lat);
grid_FSST = F_SST(lon, lat);
grid_FSSS = F_SSS(lon, lat);

if isscalar(lon)
    grid_lon = floor(lon / gap) + 1;
    grid_lat = floor((lat + 90) / gap) + 1;
    grid_mean_time = corr_mean_time(grid_lon, grid_lat);
else
    grid_lon = floor(squeeze(lon(:, 1)) / gap) + 1;
    grid_lat = floor((squeeze(lat(1, :)) + 90) / gap) + 1;
    grid_mean_time = corr_mean_time(grid_lon, grid_lat);
end

tem_distime = (datenum(messure_date, 'yyyy-mm-dd') - grid_mean_time) / 365;

DIC = grid_F0 + sst .* grid_FSST + sss .* grid_FSSS...
            + corr_b(2) .* tem_distime;

end
