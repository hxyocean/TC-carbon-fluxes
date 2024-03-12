%% This function is used to calculate the DIC at sea surface based on SST and SSS
% the recommended application scope of depth is 0-5m, but can still have 
% relative high precision with depth less than 30m.
% This program is recommended for SSS less than 30 psu
function [DIC] = fun_LIDR(SST, SSS, lon, lat, year, month, day)
% Unit:
% SST: ℃; SSS: psu; DIC: umol/kg (10-6 mol/kg)

    F = load('interpolate_cor.mat','F0','F_SST','F_SSS');
    modify_data = load('modify_data.mat');
    gap = 5.0; % the resolution of origin data map

    lon = mod(lon, 360);

    Interpolate_F0  = F.F0(lon, lat);
    Interpolate_FSST = F.F_SST(lon, lat);
    Interpolate_FSSS = F.F_SSS(lon, lat);
    Interpolate_mean_time = modify_data.mean_time(floor(lon / gap) + 1, floor((lat + 90) / gap) + 1);

    num_of_date = datenum(sprintf('%04d-%02d-%02d', year, month, day), 'yyyy-mm-dd'); % days since 0000-01-01

    DIC = Interpolate_F0 + SST * Interpolate_FSST + SSS * Interpolate_FSSS + ...
        modify_data.b(2) * (num_of_date - Interpolate_mean_time) / 365;

    
end