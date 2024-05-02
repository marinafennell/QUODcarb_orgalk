
% driver to go with QUODcarb and ORGANIG ALKALINITY stuff

load data.mat; % NEW as of Nov.11
[in] = data;
nD = length(in);

[sobs] = in(1,:); % PSU
[tobs] = in(2,:); % degC
[pobs] = in(3,:); % dbar

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
opt.pKalpha     = 0;
opt.pKbeta      = 0;

% read in GOMECC data and put into obs structure
for i = 1:nD
    obs(i).TAlpha           = 1;
    obs(i).eTAlpha          = 0.2;
    obs(i).tp(1).pKalpha    = 4.1;
    obs(i).tp(1).epKalpha   = 0.5;

    obs(i).TBeta           = 1;
    obs(i).eTBeta          = 0.2;
    obs(i).tp(1).pKbeta    = 7;
    obs(i).tp(1).epKbeta   = 0.5;

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
    % obs(i).tp(1).pKalpha = 4.38;
    % obs(i).tp(1).epKalpha = 0.5;

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
    obs(i).tp(2).T     = 20 ; %degC
    obs(i).tp(2).eT    = 0.03 ; % from cruise report
    obs(i).tp(2).P     = 0.0 ; % dbar (surface pressure for pco2)
    obs(i).tp(2).eP    = 0.07 ;
    obs(i).tp(2).pco2  = in(10,i); % (µatm)
    obs(i).tp(2).epco2 = in(10,i)*0.0021; % 0.21% relative std error (avg)
end

obs_backup = obs;
[est,obs,sys,iflag] = QUODcarb(obs,opt);

for i = 1:nD
    dTC(i) = obs(i).TC - est(i).TC;
    zscore(i,1) = dTC(i)/2.01;

    dTA(i) = obs(i).TA - est(i).TA;
    zscore(i,2) = dTA(i)/1.78;

    dph(i) = obs(i).tp(1).ph - est(i).tp(1).ph;
    zscore(i,3) = dph(i)/0.0004;

    dpco2(i) = obs(i).tp(2).pco2 - est(i).tp(2).pco2;
    zscore(i,4) = dpco2(i)/obs(i).tp(2).epco2;

    dco3(i) = obs(i).tp(1).co3 - est(i).tp(1).co3;
    zscore(i,5) = dco3(i)/obs(i).tp(1).eco3;

    % Talpha(i) = est(i).TAlpha;
    % pKalpha(i) = est(i).tp(1).pKalpha;
    % alpha(i) = est(i).tp(1).alpha;
    % halpha(i) = est(i).tp(1).halpha;
    % 
    % Tbeta(i) = est(i).TBeta;
    % pKbeta(i) = est(i).tp(1).pKbeta;
    % beta(i) = est(i).tp(1).beta;
    % hbeta(i) = est(i).tp(1).hbeta;
end

fprintf('\n')
fprintf('med zTC = %f, ', median(zscore(:,1)))
fprintf('med zTA = %f, ', median(zscore(:,2)))
fprintf('med zpH = %f, ', median(zscore(:,3)))
fprintf('\n')
fprintf('med zpCO2 = %f, ', median(zscore(:,4)))
fprintf('med zCO3 = %f ', median(zscore(:,5)))
fprintf('\n')
% fprintf('Talpha = %f, ', median(Talpha))
% fprintf('pKalpha = %f, ', median(pKalpha))
% fprintf('\n')
% fprintf('alpha = %f, ', median(alpha))
% fprintf('halpha = %f, ', median(halpha))
% fprintf('\n')
% fprintf('Tbeta = %f, ', median(Tbeta))
% fprintf('pKbeta = %f, ', median(pKbeta))
% fprintf('\n')
% fprintf('beta = %f, ', median(beta))
% fprintf('hbeta = %f, ', median(hbeta))
% fprintf('\n')


% Q5: All five input
% CT AT pH pCO2 CO3 (Q5) (fid5)
% [est,obs,sys,iflag] = QUODcarb(obs,opt);

% % TC pCO2 (Q2)(fid03)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan;         obs(i).eTA = nan;
%     obs(i).tp(2).ph = nan;   obs(i).tp(2).eph = nan;
%     obs(i).tp(2).co3 = nan;  obs(i).tp(2).eco3 = nan;
% end
% [est,obs,~,~] = QUODcarb(obs,opt);
% est03   = est;
% tp      = 3;
% fid03   = 'compare_TC_pco2.csv';
% [A]     = compare(obs,est,opt,tp,7,fid03);


