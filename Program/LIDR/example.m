%% Exmaple of LIDR method to estimate the surface DIC based on sea surface temperature and salinity;
%% Note that the program cannot automatically identify land or ocean, which should be distinguished before running.

clear;
tic;

%% single point;
sst = 30; % sea surface temperature (Celsius)
sss = 35; % sea surface salinity (psu)
lon = 180; % longitude (degree)
lat = 0;   % latitude (degree)
messure_data = '2000-01-01'; % 'yyyy-mm-dd'

DIC = fun_LIDR(sst, sss, lon, lat, messure_data);

%% mesh grid
sst = ones(360, 180) * 30;
sss = ones(360, 180) * 35;
lat = -89.5:1:89.5;
lon = 0.5:1:359.5;
[lat, lon] = meshgrid(lat, lon);
messure_data = '2000-01-01';

DIC_grid = fun_LIDR(sst, sss, lon, lat, messure_data);

%%
toc;
