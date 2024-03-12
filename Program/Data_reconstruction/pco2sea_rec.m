%% this program is use to calculate the obsvered climatoloy SSS from hycom dataset
% climatoloy SST from AVHRR dataset
tic
clear;
addpath(genpath('../jonathansharp-CO2-System-Extd-40e49f9'));
addpath(genpath('../lom310087-sup-0001-suppinfo01'));
addpath(genpath('../function'));
addpath(genpath('../LIRA2'));
%% initiation



% climatology data path
save_path = 'I:/ncs/data/data_rec/pco2sea_avhrr_hycom_smooth/';
save_path_sss = 'I:/ncs/data/data_rec/hycom_sss_025/';

%hycom data
hycom_path='H:/Hycom/';
hycom_lon_length=4500;
hycom_lon=ncread('H:\Hycom\2018\2018-01-01T000000.nc','lon');
hycom_lat=ncread('H:\Hycom\2018\2018-01-01T000000.nc','lat');
[grid_hycom_lat,grid_hycom_lon]=meshgrid(hycom_lat,hycom_lon);

hycom_lat=ncread('H:\Hycom\2021\2021-01-01T000000.nc','lat');
[grid_hycom_lat_lar,grid_hycom_lon_lar]=meshgrid(hycom_lat,hycom_lon);

hycom_lat=ncread('H:\Hycom\1993\1993-01-01T000000.nc','lat');
[grid_hycom_lat_les,grid_hycom_lon_les]=meshgrid(hycom_lat,hycom_lon);

% sst-obs data
path_ssto = 'H:/AVHRR/AVHRR.v2.1/';
lon_sst = ncread('H:/AVHRR/AVHRR.v2.1/2000/oisst-avhrr-v02r01.20000101.nc', 'lon');
lat_sst = ncread('H:/AVHRR/AVHRR.v2.1/2000/oisst-avhrr-v02r01.20000101.nc', 'lat');
lon_sst = cast(lon_sst, 'double');
lat_sst = cast(lat_sst, 'double');
[grid_lat_sst, grid_lon_sst] = meshgrid(lat_sst, lon_sst);


%-------------------------------------------------------------------------%
% DIC statistic

load('../../data/data_statistics/analysis/methods/interpolate_cor.mat','F0','F_SST','F_SSS');
modify_DIC=load('../../data/data_statistics/analysis/methods/modify.mat');
gap=5.0;

grid_lon = lon_sst;
grid_lat = lat_sst;
grid_lon(grid_lon < 0) = grid_lon(grid_lon < 0) + 360;
grid_lon = floor(grid_lon / gap) + 1;
grid_lat = floor((grid_lat + 90) / gap) + 1;
grid_lat(grid_lat > 180 / gap) = grid_lat(grid_lat > 180 / gap) - 1;
grid_lon(grid_lon > 360 / gap) = grid_lon(grid_lon > 360 / gap) - 1;
grid_mean_time = modify_DIC.mean_time(grid_lon, grid_lat);

grid_F0 = F0(grid_lon_sst, grid_lat_sst);
grid_FSST = F_SST(grid_lon_sst, grid_lat_sst);
grid_FSSS = F_SSS(grid_lon_sst, grid_lat_sst);
parsave('../../data/data_statistics/analysis/methods/grid_F.mat', ...
    grid_F0, grid_FSST, grid_FSSS, grid_mean_time);
%------------------------------------------------------------------------------------%
methods="tropical";
if methods=="tropical"
    % for tropical area
    coefficient_temp=0.0423;     coefficient_TA=-7.4;
    coefficient_DIC=8;           coefficient_sal=0.93;
elseif methods=="polar"
    % for polar area
    coefficient_temp=0.0423;     coefficient_TA=-13.3;
    coefficient_DIC=14;           coefficient_sal=1.0;
else
    % for global average
    coefficient_temp=0.0423;     coefficient_TA=-9.4;
    coefficient_DIC=10;           coefficient_sal=0.94;
end

index_tropical = (grid_lat_sst <= 40 & grid_lat_sst >= -40);
index_polor    = (grid_lat_sst >= 60 | grid_lat_sst <= -60);
index_temperate =  (grid_lat_sst > -60 & grid_lat_sst < -40) | (grid_lat_sst > 40 & grid_lat_sst < 60);

coefficient_temp = 0.0423 * ones(size(grid_lon_sst));
coefficient_DIC  = zeros(size(grid_lon_sst));
coefficient_TA  = zeros(size(grid_lon_sst));
coefficient_sal  = zeros(size(grid_lon_sst));

coefficient_DIC(index_tropical) = 8;
coefficient_DIC(index_polor)    = 14;
coefficient_DIC(index_temperate) = (abs(grid_lat_sst(index_temperate)) - 40) * (14 - 8) / 20 + 8;

coefficient_TA(index_tropical) = -7.4;
coefficient_TA(index_polor)    = -13.3;
coefficient_TA(index_temperate) = (abs(grid_lat_sst(index_temperate)) - 40) * (-13.3 + 7.4) / 20 - 7.4;

coefficient_sal(index_tropical) = 0.93;
coefficient_sal(index_polor)    = 1.0;
coefficient_sal(index_temperate) = (abs(grid_lat_sst(index_temperate)) - 40) * (1 - 0.93) / 20 + 0.93;
%------------------------------------------------------------------------------------%
% climatology pco2sea
path_pco2sea='../../data/data_statistics/pCO2 neural network/MPI-SOM_FFN_v2021_NCEI_OCADS.nc';
climate_pco2_raw = ncread(path_pco2sea,'spco2_raw');
climate_pco2_smooth = ncread(path_pco2sea,'spco2_smoothed');
climate_pco2_raw(climate_pco2_raw > 1e10) = nan;
climate_pco2_smooth(climate_pco2_smooth > 1e10) = nan;

climate_lat = ncread(path_pco2sea,'lat');
climate_lon = ncread(path_pco2sea,'lon');
climate_lon(climate_lon<0) = climate_lon(climate_lon<0)+360;
tem_climate_lon=zeros(length(climate_lon),1);
tem_climate_lon(181:360)=climate_lon(1:180);
tem_climate_lon(1:180)=climate_lon(181:360);
climate_lon(:)=tem_climate_lon(:);

climate_lon = [climate_lon(end)-360; climate_lon; climate_lon(1)+ 360];
[grid_lat_climate, grid_lon_climate] = meshgrid(climate_lat, climate_lon);


tem_data=zeros(size(climate_pco2_raw));
tem_data(1:180,:,:)=climate_pco2_raw(181:360,:,:);
tem_data(181:360,:,:)=climate_pco2_raw(1:180,:,:);
climate_pco2_raw(:)=tem_data(:);


tem_data=zeros(size(climate_pco2_smooth));
tem_data(1:180,:,:)=climate_pco2_smooth(181:360,:,:);
tem_data(181:360,:,:)=climate_pco2_smooth(1:180,:,:);
climate_pco2_smooth(:)=tem_data(:);


climate_begintime='1982-01'; % the begin time of climate

clear tem_climate_lon tem_data
%% the climatology SST and SSS

time_begin_sss='1992-12';
time_end_sss='2020-12';
time_begin_ssto='1992-12';
time_end_ssto='2020-12';

time_begin = datestr(min(datenum(time_begin_ssto), datenum(time_begin_sss)), 'yyyy-mm');
time_end   = datestr(max(datenum(time_end_ssto), datenum(time_end_sss)), 'yyyy-mm');
num_month=(str2double(time_end(1:4))-str2double(time_begin(1:4)))*12+...
    (str2double(time_end(6:7))-str2double(time_begin(6:7)))+1;

delete(gcp('nocreate'));         % stop the process before start a new run.
numCore = feature('numcores');   % get the maxmium core num of PC.
parpool(numCore-7);                % start parpool.

parfor i = 1:num_month
    current_month = mod(str2double(time_begin(6:7)) + i - 1 - 1, 12) + 1; 
    next_month = mod(str2double(time_begin(6:7)) + i - 1 - 1 + 1, 12) + 1; 
    current_year  = str2double(time_begin(1:4)) + ...
        floor((str2double(time_begin(6:7)) + i - 1 - 1) / 12);
    next_year = str2double(time_begin(1:4)) + ...
        floor((str2double(time_begin(6:7)) + i - 1 - 1 + 1) / 12);
    current_month = sprintf('%02d', current_month);
    current_year  = sprintf('%04d', current_year);
    next_month = sprintf('%02d', next_month);
    next_year  = sprintf('%04d', next_year);
    current_time  = [current_year, '-', current_month];
    next_time  = [next_year, '-', next_month];
    tem_num  = datenum([next_time,'-01']) - datenum([current_time,'-01']);

    tem_index_climate=(str2double(current_year)-str2double(climate_begintime(1:4))) * 12 +...
            (str2double(current_month)-str2double(climate_begintime(6:7))) + 1;

    climatology_pco2sea = interp2(grid_lat_climate, grid_lon_climate, ...
        [climate_pco2_smooth(end,:,tem_index_climate);...
        climate_pco2_smooth(:,:,tem_index_climate); ...
        climate_pco2_smooth(1,:,tem_index_climate)],...
        grid_lat_sst, grid_lon_sst, 'linear');

    if exist([save_path, current_time(1:4), '/', current_time, '-01_pco2sea.mat'], "file")
        continue;
    end
    disp(['calculating ',current_time])
%     climate_index=(str2double(current_year)-str2double(climate_begintime(1:4))) * 12 +...
%         (str2double(current_month)-str2double(climate_begintime(6:7)))+1;
%     tem_climate_pco2_raw=climate_pco2_raw(:,:,climate_index);
%     tem_climate_pco2_smooth=climate_pco2_smooth(:,:,climate_index);
    

    SST     = nan(length(lon_sst), length(lat_sst), tem_num,'single');
    SSS     = nan(length(lon_sst), length(lat_sst), tem_num,'single');
    DIC     = nan(length(lon_sst), length(lat_sst), tem_num,'single');
    TA     = nan(length(lon_sst), length(lat_sst), tem_num,'single');
    

    for j = 1:tem_num
        
        tem_time_current = datestr(datenum([current_time,'-01']) + j - 1, ...
            'yyyy-mm-dd');

        %=======================================================================================%
        % hycom sss
        tem_hycom_name_00=[hycom_path,tem_time_current(1:4),'/',tem_time_current,'T000000.nc'];
        if ~exist(tem_hycom_name_00,'file')
            tem_str = datestr(datenum([current_time,'-01']) + j - 2, 'yyyy-mm-dd');
            tem_hycom_name_00 = [hycom_path,tem_str(1:4),'/',tem_str,'T210000.nc'];
        end
        tem_hycom_name_12=[hycom_path,tem_time_current(1:4),'/',tem_time_current,'T120000.nc'];
        if ~exist(tem_hycom_name_12,'file')
            tem_hycom_name_12=[hycom_path,tem_time_current(1:4),'/',tem_time_current,'T090000.nc'];
            if ~exist(tem_hycom_name_12,'file')
                tem_hycom_name_12=[hycom_path,tem_time_current(1:4),'/',tem_time_current,'T150000.nc'];
            end
        end
        if ~exist(tem_hycom_name_00,'file')
            tem_hycom_name_00 = tem_hycom_name_12;
        end
        if ~exist(tem_hycom_name_00,'file') && ~exist(tem_hycom_name_12,'file')
            disp([current_time,' missing']);
            continue;
        end
        tem_hycom_lon_00=ncread(tem_hycom_name_00,'lon');
        tem_hycom_lon_12=ncread(tem_hycom_name_12,'lon');
        
        
        SSS_00=ncread(tem_hycom_name_00,'salinity');        
        SSS_12=ncread(tem_hycom_name_12,'salinity');

        if tem_hycom_lon_00(1) < 0
            tem_data=SSS_00;
            tem_data(1:2250,:)=SSS_00(2251:end,:);
            tem_data(2251:end,:)=SSS_00(1:2250,:);
            SSS_00=tem_data;
        end
        if tem_hycom_lon_12(1) < 0
            tem_data=SSS_12;
            tem_data(1:2250,:)=SSS_12(2251:end,:);
            tem_data(2251:end,:)=SSS_12(1:2250,:);
            SSS_12=tem_data;
        end
        sss_025 = nanmean(cat(3,SSS_00, SSS_12), 3);
        
        if datenum([current_time,'-01']) + j - 1 >= datenum('2020-02-01')
            sss_025 = interp2(grid_hycom_lat_lar, grid_hycom_lon_lar, sss_025, ...
                grid_lat_sst, grid_lon_sst, "linear");
        elseif datenum([current_time,'-01']) + j - 1 <= datenum('1993-12-31')
            sss_025 = interp2(grid_hycom_lat_les, grid_hycom_lon_les, sss_025, ...
                grid_lat_sst, grid_lon_sst, "linear");
        else
            sss_025 = interp2(grid_hycom_lat, grid_hycom_lon, sss_025, ...
                grid_lat_sst, grid_lon_sst, "linear");
        end
        SSS(:, :, j) = sss_025;
        
        %=======================================================================================%
        % AVHRR sst
        tem_file_name = [path_ssto, tem_time_current(1:4), '/', ...
            'oisst-avhrr-v02r01.', tem_time_current([1,2,3,4,6,7,9,10]),'.nc'];
        sst_025 = ncread(tem_file_name, 'sst');
        SST(:, :, j) = sst_025;

        %=======================================================================================%
        % DIC
        
        tem_distime = (datenum(tem_time_current) - grid_mean_time)/365;
        
        DIC_025 = grid_F0 + sst_025 .* grid_FSST + sss_025 .* grid_FSSS...
            + modify_DIC.b(2) .* tem_distime;
        
        DIC(:, :, j) = DIC_025;
        %=======================================================================================%
        % TA
        [TA_025, uncret] = LIAR2([grid_lon_sst(:), grid_lat_sst(:),...
            zeros(length(grid_lat_sst(:)),1)],...
            [sss_025(:), sst_025(:)], [1, 2], 'Equations', 8);
        TA_025 = reshape(TA_025,size(sst_025));

        TA(:, :, j) = TA_025;

    end

    %tem_SSS(tem_SSS<30)=NaN;
    %=====================================================================================%
    sst_025 = mean(SST, 3);
    sss_025 = mean(SSS, 3);
    DIC_025 = mean(DIC, 3);
    TA_025 = mean(TA, 3);

    for j = 1:tem_num
        tem_time_current = datestr(datenum([current_time,'-01']) + j - 1, ...
            'yyyy-mm-dd');
        pco2_sea = climatology_pco2sea + climatology_pco2sea .* ...
                    (coefficient_temp.*(SST(:, :, j) - sst_025) ... % temperature
                    +coefficient_DIC .* (DIC(:, :, j) - DIC_025) ./ DIC_025 ... % DIC
                    +coefficient_TA .* (TA(:, :, j) - TA_025) ./ TA_025 ... % TA
                    +coefficient_sal .* (SSS(:, :, j) - sss_025) ./ sss_025); % SSS;
        if ~exist([save_path, tem_time_current(1:4)], "dir")
            mkdir([save_path, tem_time_current(1:4)]);
        end
        parsave([save_path, tem_time_current(1:4), '/', tem_time_current([1,2,3,4,6,7,9,10]), '_pco2sea.mat'], pco2_sea);
        
        sss = SSS(:, :, j);
        if ~exist([save_path_sss, tem_time_current(1:4)], "dir")
            mkdir([save_path_sss, tem_time_current(1:4)]);
        end
        parsave([save_path_sss, tem_time_current(1:4), '/', tem_time_current([1,2,3,4,6,7,9,10]), '_sss.mat'], sss);

    end
    disp([current_time,' completed!'])
end

toc;
