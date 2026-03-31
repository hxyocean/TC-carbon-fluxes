
%% This program is used to examine the influence of the nonlinearity of pCO2

%%
tic;
clear;

addpath(genpath('../jonathansharp-CO2-System-Extd-40e49f9')); % import the function of CO2SYS


%% the parameter of CO2SYS
SCALE  = 1; % Total pH scale
K1K2   = 10; % Leuker et al (2000) K1K2
SO4    = 1; % Dickson (1990) KSO4
KF     = 2; % Perez & Fraga (1987) KF
BOR    = 2; % Lee et al (2010) TB
SI = 3.9;
PSO4 = 0.59;

%% parameter for PRT vertical mixing model
Ri = 1;
a = 0.0002;
g = 9.81;
f = 5e-5;
Cd = 1.5e-3;
V = 40;
rho_sea = 1.025e3;
rho_air = 1.25;
t = rho_air * Cd * V * V;

%%
% h_pre: pre-storm mixed layer depth;
% h_aft: post-storm mixed layer depth

h_pre = 50; 

% the F_ denote the vertical gradient of corresponding parameters below
% the mixed layer; v_pre_ denote the pre-storm surface value of 
% corresponding parameters 

F_DIC = -0.4;
F_TA  = -0.25;
F_SST = 0.05;
F_SSS = -0.0;


v_pre_DIC = 2000;
v_pre_TA  = 2300;
v_pre_SST = 30;
v_pre_SSS = 35;


num = 101;
var_DIC = linspace(-200, 200, num);
var_TA  = linspace(-200, 200, num);
var_SST = linspace(-5, 5, num);
var_SSS = linspace(-5, 5, num);
var_F_DIC = linspace(-0.1, 0.1, num);
var_F_TA = linspace(-0.1, 0.1, num);
var_F_SST = linspace(-0.01, 0.01, num);
var_F_SSS = linspace(-0.005, 0.005, num);



dpco2_all = [];


for m = 1:8
    dpco2_linear = zeros(num, 1);
    dpco2_nolinear = zeros(num, 1);

    for n = 1:num
        if m == 1
            tem_v_pre_DIC = v_pre_DIC + var_DIC(n);
            tem_v_pre_TA = v_pre_TA;
            tem_v_pre_SST = v_pre_SST;
            tem_v_pre_SSS = v_pre_SSS;
            tem_F_DIC = F_DIC;
            tem_F_TA = F_TA;
            tem_F_SST = F_SST;
            tem_F_SSS = F_SSS;
        elseif m == 2
            tem_v_pre_DIC = v_pre_DIC;
            tem_v_pre_TA = v_pre_TA + var_TA(n);
            tem_v_pre_SST = v_pre_SST;
            tem_v_pre_SSS = v_pre_SSS;
            tem_F_DIC = F_DIC;
            tem_F_TA = F_TA;
            tem_F_SST = F_SST;
            tem_F_SSS = F_SSS;            
        elseif m == 3
            tem_v_pre_DIC = v_pre_DIC;
            tem_v_pre_TA = v_pre_TA;
            tem_v_pre_SST = v_pre_SST + var_SST(n);
            tem_v_pre_SSS = v_pre_SSS;   
            tem_F_DIC = F_DIC;
            tem_F_TA = F_TA;
            tem_F_SST = F_SST;
            tem_F_SSS = F_SSS;  
        elseif m == 4
            tem_v_pre_DIC = v_pre_DIC;
            tem_v_pre_TA = v_pre_TA;
            tem_v_pre_SST = v_pre_SST;
            tem_v_pre_SSS = v_pre_SSS + var_SSS(n);   
            tem_F_DIC = F_DIC;
            tem_F_TA = F_TA;
            tem_F_SST = F_SST;
            tem_F_SSS = F_SSS;  
        elseif m == 5
            tem_v_pre_DIC = v_pre_DIC;
            tem_v_pre_TA = v_pre_TA;
            tem_v_pre_SST = v_pre_SST;
            tem_v_pre_SSS = v_pre_SSS;   
            tem_F_DIC = F_DIC + var_F_DIC(n);
            tem_F_TA = F_TA;
            tem_F_SST = F_SST;
            tem_F_SSS = F_SSS;  
        elseif m == 6
            tem_v_pre_DIC = v_pre_DIC;
            tem_v_pre_TA = v_pre_TA;
            tem_v_pre_SST = v_pre_SST;
            tem_v_pre_SSS = v_pre_SSS;   
            tem_F_DIC = F_DIC;
            tem_F_TA = F_TA + var_F_TA(n);
            tem_F_SST = F_SST;
            tem_F_SSS = F_SSS;  
        elseif m == 7
            tem_v_pre_DIC = v_pre_DIC;
            tem_v_pre_TA = v_pre_TA;
            tem_v_pre_SST = v_pre_SST;
            tem_v_pre_SSS = v_pre_SSS;   
            tem_F_DIC = F_DIC;
            tem_F_TA = F_TA;
            tem_F_SST = F_SST + var_F_SST(n);
            tem_F_SSS = F_SSS;  
        else 
            tem_v_pre_DIC = v_pre_DIC;
            tem_v_pre_TA = v_pre_TA;
            tem_v_pre_SST = v_pre_SST;
            tem_v_pre_SSS = v_pre_SSS;   
            tem_F_DIC = F_DIC;
            tem_F_TA = F_TA;
            tem_F_SST = F_SST;
            tem_F_SSS = F_SSS + var_F_SSS(n);  
        end

        

        N2 = g * a * F_SST;
        h_aft = h_pre * h_pre + sqrt(power(h_pre, 4) + 32 * Ri * t * t / N2 / rho_sea / rho_sea / f / f);
        h_aft = sqrt(h_aft / 2);
        depth = 0:0.1:h_aft;

        %%
        pro_DIC = zeros(size(depth));
        pro_TA = zeros(size(depth));
        pro_SST = zeros(size(depth));
        pro_SSS = zeros(size(depth));
        pro_pco2 = zeros(size(depth));


        for i = 1:length(depth)
            if depth(i) <= h_pre
                
                pro_SST(i) = tem_v_pre_SST;
                pro_SSS(i) = tem_v_pre_SSS;
                pro_DIC(i) = tem_v_pre_DIC;
                pro_TA(i) = tem_v_pre_TA;
            else
                pro_SST(i) = tem_v_pre_SST - tem_F_SST * (depth(i) - h_pre);
                pro_SSS(i) = tem_v_pre_SSS - tem_F_SSS * (depth(i) - h_pre);
                pro_DIC(i) = tem_v_pre_DIC - tem_F_DIC * (depth(i) - h_pre);
                pro_TA(i) = tem_v_pre_TA - tem_F_TA * (depth(i) - h_pre);
            end
            
        end
        
        [C,Headers,Niceheaders] = CO2SYS(pro_TA,pro_DIC,1,2,pro_SSS,pro_SST, ...
                    nan,0,nan,SI,PSO4,0,0,SCALE,K1K2,SO4,KF,BOR);
        pro_pco2 = C(:, 3);


        avg_DIC = 0;
        avg_TA = 0;
        avg_SST = 0;
        avg_SSS = 0;
        avg_pco2 = 0;
        
        for i = 1:length(depth)-1
            avg_DIC = avg_DIC + (pro_DIC(i+1) + pro_DIC(i)) / 2 * (depth(i+1) - depth(i));
            avg_TA = avg_TA + (pro_TA(i+1) + pro_TA(i)) / 2 * (depth(i+1) - depth(i));
            avg_SST = avg_SST + (pro_SST(i+1) + pro_SST(i)) / 2 * (depth(i+1) - depth(i));
            avg_SSS = avg_SSS + (pro_SSS(i+1) + pro_SSS(i)) / 2 * (depth(i+1) - depth(i));
            avg_pco2 = avg_pco2 + (pro_pco2(i+1) + pro_pco2(i)) / 2 * (depth(i+1) - depth(i));
        
            
        end
        
        avg_DIC = avg_DIC / (depth(end) - depth(1));
        avg_TA = avg_TA / (depth(end) - depth(1));
        avg_SST = avg_SST / (depth(end) - depth(1));
        avg_SSS = avg_SSS / (depth(end) - depth(1));
        avg_pco2 = avg_pco2 / (depth(end) - depth(1));
        
        [C,Headers,Niceheaders] = CO2SYS(avg_TA,avg_DIC,1,2,avg_SSS,avg_SST, ...
                    nan,0,nan,SI,PSO4,0,0,SCALE,K1K2,SO4,KF,BOR);
        
        avg_pco2_nolinear = C(3);

        dpco2_linear(n) = avg_pco2 - pro_pco2(1);
        dpco2_nolinear(n) = avg_pco2_nolinear - pro_pco2(1);


    end

    dpco2_all = [dpco2_all, [dpco2_linear, dpco2_nolinear]];
end


%%

save('../intermediate_file/pco2.mat', "dpco2_all", ...
    "var_F_SST", "var_F_TA", "var_F_DIC", "var_F_SSS", ...
    "var_SST", "var_TA", "var_DIC", "var_SSS", ...
    "v_pre_SST", "v_pre_DIC", "v_pre_TA", "v_pre_SSS", ...
    "F_SST", "F_DIC", "F_TA", "F_SSS");



toc;
