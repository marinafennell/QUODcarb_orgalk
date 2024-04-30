
% new driver for opt.turnoff.pK1 testing with just tp(1)

load data.mat;
[in] = data;
nD = 5; % length(in);

% choose options for opt structure
opt.K1K2 = 16; % option for K1K2 formulation
opt.KSO4 = 1;  % option for KSO4 formulation
opt.KF   = 2;  % option for KF formulation
opt.TB   = 2;  % option for TB formulation
opt.phscale  = 1;  % 1 = tot, 2 = sws, 3 = free, 4 = NBS
opt.printcsv = 0;  % print est to CSV? 1 = on , 0 = off
% opt.fname    = 'QUODcarb_output.csv'; % don't need it if printcsv is off
opt.fname    = 'Q5_pK1_mar14.csv';
opt.co2press = 1; % 1 = on, 0 = off
opt.Revelle  = 0; % 1 = on, 0 = off 
opt.printmes = 0; % 1 = on, 0 = off

opt.turnoff.TB  = 0; % 1 = on (no TB formulation used)
opt.turnoff.pK1 = 0;
opt.pKalpha = 1;

% load output_mat_files/K16/est09; % calculated pCO2 at 25C
% for FP
obs(1).tp(1).pco2 = 363.1889;
obs(1).tp(1).epco2 = 0.0021*obs(1).tp(1).pco2;
obs(2).tp(1).pco2 = 356.9050;
obs(2).tp(1).epco2 = 0.0021*obs(2).tp(1).pco2;
obs(3).tp(1).pco2 = 947.6132;
obs(3).tp(1).epco2 = 0.0021*obs(3).tp(1).pco2;
obs(4).tp(1).pco2 = 955.3171;
obs(4).tp(1).epco2 = 0.0021*obs(4).tp(1).pco2;
obs(5).tp(1).pco2 = 961.0366;
obs(5).tp(1).epco2 = 0.0021*obs(5).tp(1).pco2;

for i = 1:nD
    
    % measurements that are independent of (T,P)
    obs(i).TC    = in(5,i); % (umol/kg)
    obs(i).eTC   = 2.01;    % TC error ±2.01 umol/kg
    obs(i).TA    = in(6,i);
    obs(i).eTA   = 1.78; % TA error ±2.01 umol/kg
    obs(i).sal   = in(1,i); % PSU
    obs(i).esal  = 0.001; % 1 new of 1/23 old = 0.002
    % nutrients P and Si also independent of (T,P)
    obs(i).TP    = in(7,i);
    obs(i).eTP   = in(7,i)*0.003; % 0.30% meas precision NEW 4/17/24
    obs(i).TSi   = in(8,i);
    obs(i).eTSi  = in(8,i)*0.0031; % 0.31% meas uncertainty NEW 4/17/24

    % second(T,P)-dependent measurement
    obs(i).tp(1).T    = 25 ; % degC
    obs(i).tp(1).eT   = 0.05 ; % from cruise report
    obs(i).tp(1).P    = 0.0 ; %in(i+ad,1); % NOT in situ
    obs(i).tp(1).eP   = 0.07 ;
    obs(i).tp(1).ph   = in(9,i); % total scale
    obs(i).tp(1).eph  = 0.0004 ;
    obs(i).tp(1).co3  = in(11,i); % (µmol/kg)
    obs(i).tp(1).eco3 = in(11,i)*0.02; % 2% from Jon Sharp NEW 1/25/24

    % obs(i).tp(1).pco2 = est09(i).tp(2).pco2;
    % obs(i).tp(1).epco2 = est09(i).tp(2).epco2;
end

[est,obs,sys,iflag] = QUODcarb(obs,opt);





