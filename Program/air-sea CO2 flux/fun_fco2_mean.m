% this function is used to calculate the air-sea flux of co2
%  based on the wind speed and air-sea pCO2 difference
%  temperature is needed to determine the Sc (k = 0.251*u^2*(Sc/660)^-0.5)
% the equation is FCO2 = k * KH * ΔpCO2
% the units of dpco2 is uatm 
function [fco2, k, k_std, kh] = fun_fco2_mean(wind_sq, dpco2, temp, sal)
    
    % wind_sq: wind speed (m/s)
    % dpco2: air-sea pCO2 difference (uatm)

    A1 = -58.0931; %-60.2409
    A2 = 90.5069; %93.4517
    A3 = 22.2940; %23.3585
    B1 = 0.027766; %0.023517
    B2 = -0.025888; %-0.023656
    B3 = 0.0050578; %0.0047036

    %% calculate the Sc k and kh
    % Sc is the function of temperature (Wanninkhof, 2014)
    % kh is the solubility of CO2 that is a function of salinity and temperature (Weiss, 1974)
    temp(temp<-2)=nan;
    sal(sal<30)=nan;
    
    Sc = 2116.8 - 136.25 .* temp ...
        + 4.7353 .* temp .* temp...
        - 0.092307 .* temp .* temp.* temp...
        + 0.0007555 .* temp .* temp.* temp .* temp;

    temp=temp+273.15;
    kh = exp(A1 + A2.*(100 ./ temp) + A3.*log(temp./100) ...
        + sal .* (B1 + B2 .* (temp./100) + B3 .* (temp./100) .* (temp./100)));%units mol/(L*atm)


    k1 = (0.233 * wind_sq.*wind_sq + 0.333 .* wind_sq).*((Sc/600).^(-0.5)); %units cm/h
    k2 = 0.251.*wind_sq.*wind_sq.*((Sc/660).^(-0.5)); %units cm/h
    k3 = 0.0283.*wind_sq.*wind_sq.*wind_sq.*((Sc/660).^(-0.5)); %units cm/h

    k = (k1+k2+k3)/3;
    k_std = sqrt(((k1 - k) .* (k1 - k) + (k2 - k) .* (k2 - k) + (k3 - k) .* (k3 - k)) / 3) / sqrt(3);
    
    %% compute the fco2
    
    fco2 = k.*kh.*dpco2 * 1e-6; % units:  mol/m2/hour * 10
    fco2 = fco2 * 10 * 24;  % units:  mol/m2/day

end