
% driver to go with QUODcarb and ORGANIG ALKALINITY stuff

load data.mat; % NEW as of Nov.11
[in] = data;
nD = 5; %length(in);

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
opt.pKalpha     = 1;    
opt.pKbeta      = 0;

% read in GOMECC data and put into obs structure
for i = 1:nD
    % obs(i).TAlpha           = 1;
    % obs(i).eTAlpha          = 0.05;
    % obs(i).tp(1).pKalpha    = 4.0; % epK is 10% on K

    % obs(i).TBeta            = 10;
    % obs(i).eTBeta           = 0.5;
    % obs(i).tp(1).pKbeta     = 7.8; % epK is 10% on K, built in

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

    % first (T,P)-dependent measurement
    % obs(i).tp(1).T  = in(2,i); % deg C, CTD temp
    % obs(i).tp(1).eT = 0.02; % ±0.02 degC
    % obs(i).tp(1).P  = in(3,i); % dbar
    % obs(i).tp(1).eP = 0.63; % (max) ± 0.63 dbar

    % second(T,P)-dependent measurement
    obs(i).tp(1).T    = 25 ; % degC
    obs(i).tp(1).eT   = 0.05 ; % from cruise report
    obs(i).tp(1).P    = 0.0 ; %in(i+ad,1); % NOT in situ
    obs(i).tp(1).eP   = 0.07 ;
    obs(i).tp(1).ph   = in(9,i); % total scale
    obs(i).tp(1).eph  = 0.0004 ;
    obs(i).tp(1).co3  = in(11,i); % (µmol/kg)
    obs(i).tp(1).eco3 = in(11,i)*0.02; % 2% from Jon Sharp NEW 1/25/24

    % third (T,P)-dependent measurement
    obs(i).tp(2).T     = 20 ; % degC
    obs(i).tp(2).eT    = 0.03 ; % from cruise report
    obs(i).tp(2).P     = 0.0 ; % dbar (surface pressure for pco2)
    obs(i).tp(2).eP    = 0.07 ;
    obs(i).tp(2).pco2  = in(10,i); % (µatm)
    obs(i).tp(2).epco2 = in(10,i)*0.0021; % 0.21% relative std error (avg)
end
obs_backup = obs;
% keyboard

% [est,obs,sys,iflag] = QUODcarb(obs,opt);

Tin = [0.1; 1; 5; 10; 20; 50; 100; 200; 500; 1000; 2000; 3000];
pKin = 3.4:0.5:8.4;
[X,Y] = meshgrid(Tin,pKin);
Z = nan*X;
Zind = nan*X;

% for k = 1:length(pKin)
%     for j= 1:length(Tin)
for j = 1:length(X(:))
    obs = obs_backup;
    x = X(j);
    y = Y(j);
    for i = 1:nD
        % obs(i).TAlpha           = Tin(j);
        % obs(i).eTAlpha          = 0.05*Tin(j); % eTAlpha is 5%
        % obs(i).tp(1).pKalpha    = pKin(k); % epK is 10% on K
        obs(i).TAlpha           = x;
        obs(i).eTAlpha          = 0.05*x; % eTAlpha is 5%
        obs(i).tp(1).pKalpha    = y; % epK is 10% on K
    end
    [est,~,~,~]     = QUODcarb(obs,opt);
    fname = sprintf('output_mat_files/all132/est_%d',j);
    save(fname,'est');

    for i = 1:nD
        f(i) = est(i).f;
    end
    [ind,mn] = min(f);
    Z(j) = mn;
    Zind(j) = ind;
    clear est;
%     end
% end
end








