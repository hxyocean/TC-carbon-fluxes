%% exmaple for the calculation of FCO2 anomaly induced by TC
tic;
clear;
%% pre-storm value

wind_sq = 10;                % sea surface wind speed (m/s)
pco2_sea = 300;              % seawater pco2 at the air-sea interface (uatm)
pco2_air = 270;              % atmospheric pco2 at the air-sea interface (uatm)
dpco2 = pco2_sea - pco2_air; % the difference of pco2sea and pco2air (uatm)
temp = 30;                   % sea surface temperature (Celsius)
sal = 35;                    % sea surface salinity (psu)

[fco2_pre, k, k_std, kh] = fun_fco2_mean(wind_sq, dpco2, temp, sal);


%% TC forcing

wind_sq = 10 + 30.;   % abrupt wind increase
pco2_sea = 300 - 10.; % pco2sea decrease
pco2_air = 270 - 5.;  % pco2sea decrease due to low pressure
dpco2 = pco2_sea - pco2_air;
temp = 30 - 2.;       % sea surface temperature decrease due to cold wake
sal = 35;

[fco2_during, k, k_std, kh] = fun_fco2_mean(wind_sq, dpco2, temp, sal);


% caluclate the FCO2 anomaly induced by TC (unit: mol/m2/day)
delt_fco2 = fco2_during - fco2_pre;


toc;