function [est,obs,sys,iflag] = QUODcarb(obs,opt) % ORG alk
%
% OUTPUT:
%   est := posterior estimates of co2-system variables, equilibrium constants, and precisions
%   obs := same as input except that the pK's have been added
% iflag := 0 if solver converged to specified accuracy 
%          1 after reaching maximum number of iterations without converging
%          2 if it was one order of magnitude away from converging at
%               maximum iteration
%   CSV := optional CSV output, if turned on in options
%
% INPUT:
%   obs  := co2-system measured quantities with precisions 
%   opt  := solver options input by user 
%
% -------------------------------------------------------------------------
%
% SYNTAX example:
%   obs.sal         = salinity;     (PSU)           
%   obs.esal        = sal_error;    (±sigma)        
%   obs.TC          = total_c;      (umol/kg-SW)    
%   obs.eTC         = TC_error;     (±sigma)        
%   obs.TA          = total_alk;    (umol/kg-SW)    
%   obs.eTA         = alk_error;    (±sigma)        
%   obs.tp(1).T     = temp;         (deg C)         
%   obs.tp(1).eT    = temp_error;   (±sigma)        
%   obs.tp(1).P     = pressure;     (dbar, 0 = surface)          
%   obs.tp(1).eP    = pres_error;   (±sigma)     
%   obs.tp(1).ph    = ph_meas;      
%   obs.tp(1).eph   = ph_error;     (±sigma)
%
%   opt.K1K2        = 10;           % (Lueker et al 2000)
%   opt.KSO4        = 1;            % (Dickson et al 1990a) 
%   opt.KF          = 2;            % (Perez and Fraga 1987
%   opt.TB          = 2;            % (Lee et al. 2010)
%   opt.phscale     = 1;            % (1=tot, 2=free, 3=sws, 4=nbs)
%   opt.printcsv    = 1;            % (1=on, 0=off)
%   opt.fname       = 'output.csv'; % (CSV filename)
%   opt.printmes    = 1;            % (1=on, 0=off)
%   opt.co2press    = 1;            % (1=on, 0=off)
%   opt.Revelle     = 1;            % (1=on, 0=off)
%
%--------------------------------------------------------------------------
% 
% INPUT OPTIONS:
%   opt.K1K2  -> choice of K1 and K2 formulation
%           1 = Roy et al, 1993
%           2 = Goyet & Poisson, 1989
%           3 = Hansson, 1973          REFIT by Dickson and Millero, 1987
%           4 = Mehrbach et al., 1973  REFIT by Dickson and Millero, 1987
%           5 = Hansson, 1973 and Mehrbach, 1973 
%                                      REFIT by Dickson and Millero, 1987
%        x(6) = x(GEOSECS)            ~NOT AVAILABLE IN QUODcarb~
%        x(7) = x(Peng)               ~NOT AVAILABLE IN QUODcarb~
%        x(8) = x(Millero, 1979)      ~NOT AVAILABLE IN QUODcarb~
%           9 = Cai and Wang, 1998
%          10 = Lueker et al., 2000    (DEFAULT)
%          11 = Mojica Prieto and Millero, 2002
%          12 = Millero et al., 2000
%          13 = Millero et al., 2002
%          14 = Millero et al., 2006
%          15 = Waters, Millero, and Woosley, 2014
%          16 = Sulpis et al., 2020
%          17 = Schockman and Byrne, 2021
%
%   opt.KSO4  -> choice of KSO4 formulation
%           1 = Dickson (1990a)         (DEFAULT)
%           2 = Khoo et al., 1977
%           3 = Waters and Millero, 2013
%
%   opt.KF    -> choice of KF formulation
%           1 = Dickson and Riley, 1979
%           2 = Perez and Fraga, 1987   (DEFAULT)
%
%   opt.TB    -> choice of total borate formulation
%           1 = Uppstrom, 1979
%           2 = Lee et al., 2010        (DEFAULT)
%
%   opt.co2press -> turn on or off the pressure dependencies for K0 and
%           pCO2 to fCO2 fugacity factor (p2f)
%           1 = on                      (DEFAULT)
%           2 = off
%
%--------------------------------------------------------------------------
%
% OUTPUT:
%   est ->  'est' structure with best estimate contains:
%               1. p(value) and p(error) where p(x) = -log10(x)
%               2. value and average error about the value in 'q'
%                   where q(x) = x^(-10)
%               3. upper and lower bounds in 'q' space, not symmetric
%                   about the value in 'q' space
%   csv ->  csv file with most of est populated in a spreadsheet, 
%                 contains column headers with labels and units                 
%                    -does not include upper and lower errors
%
%--------------------------------------------------------------------------
%
% Changes? -> the only things you may want to change are: 
%               1. tolerance level of Newton solver -> line 111
%               2. Max Iteration number -> MAXIT in newtn.m
%
%--------------------------------------------------------------------------


    opt = check_opt(opt);                   % check opt structure
    sys = mksys(obs(1),opt.phscale,opt.pKalpha,opt.pKbeta);
                                            % mksys creates indexing 
    nD  = length(obs);  
    nv  = size(sys.K,2);

    % populate obs, yobs, wobs at each datapoint
    [obs,yobs,wobs,sys] = parse_input(obs,sys,opt,nD);

    for i = 1:nD % loop over the full data set

        z0      = init(opt,yobs(i,:),sys);   % initialize
        tol     = 1e-6;                      % tolerance

        % negative of the log of the posterior 
        % aka the log improbability (limp for short)
        fun = @(z) limp(z,yobs(i,:),wobs(i,:),obs(i),sys,opt);
        
        % find the maximum of the posterior probability 
        [zhat,J,iflag(i)] = newtn(z0,fun,tol);
        if (iflag(i) ~=0) && (opt.printmes ~= 0)
            fprintf('Newton''s method iflag = %i at i = %i \n',iflag(i),i);
        end
        [~,~,f] = limp(zhat,yobs(i,:),wobs(i,:),obs(i),sys,opt);
        % residual f value, to tack onto est
        % keyboard
        % calculate the marginalized posterior uncertainty using Laplace's approximation
        C = inv(J);
        C = C(1:nv,1:nv);
        sigx = sqrt(full(diag(C)));
        if (opt.printmes ~= 0)
            if (sum(isnan(sigx)) > 0) || (sum(isinf(sigx)) > 0) 
            fprintf('NaN found in output means faulty run. i = %i\n',i)
            end
        end

        % populate est
        [est(i)] = parse_output(zhat,sigx,sys,f);    

        % calculate the Revelle factor if opt.Revelle = 1
        if opt.Revelle == 1
            for j = 1:length(sys.tp(:))
                % Revelle
                ifree   = sys.tp(j).ifree;

                ei      = zeros(length(ifree),1);
                ei(1)   = 1;
                jac     = sys.tp(j).dcdx_pTAfixed(zhat(ifree));
                z       = ei - ( jac.' ) * ( ( jac * jac.' ) \ ( jac*ei ) );
                est(i).tp(j).Revelle = z(2)/z(1);

                % dpfCO2dpTA (similar to Revelle but TC held fixed)
                jfree   = sys.tp(j).jfree;
                ej      = zeros(length(jfree),1);

                ej(1)   = 1;
                jac     = sys.tp(j).dcdx_pTCfixed(zhat(jfree)) ;
                z       = ej - ( jac.' ) * ( ( jac * jac.' ) \ ( jac*ej ) );
                est(i).tp(j).dpfco2dpTA = z(2)/z(1);
            end
        end
    end

    % PrintCSV if opt.printcsv = 1 using filename opt.fname
    PrintCSV(est,obs,iflag,opt);
end

% -------------------------------------------------------------------------
% Subfunctions
% -------------------------------------------------------------------------

function [g,H,f] = limp(z,y,w,obs,sys,opt)
% [g,H,f] = limp(z,y,w,obs,sys,opt)
%
% Negative log probability for the co2 system  a.k.a. log improbability i.e., limp!
%
% INPUT:
%
%   y  := measured components, non-measured components set to NaN
%   w  := measurement precisions, same size and order as y, non-measured components set to anything, they are ignored
% gpK  := first order derivatives of pK wrt T,S,P
% ggpK := second order derivatives of pK wrt T,S,P

% OUTPUT:
%
%   f := limp
%   g := grad f w.r.t. x
%   h := hessian of f w.r.t. x

    p   = sys.p;    
    q   = sys.q;
    M   = sys.M;    
    K   = sys.K;
    nrk     = size(K,1);
    nTP     = length(sys.tp); 
    nv      = size(M,2);
    x       = z(1:nv);      % thermodynamic system state, PP*x is modeled measured vars

    % fill ypK, gypK, ggypK, and with associated calculated pK, gpK, and ggpK values
    % update the pK values based on the new estimate of (T,P)
    [y, gy, ggy ] = update_y(y,x,obs,sys,opt);
    % keyboard
    % Make a vector of measured quantities    
    id  = find(~isnan(y));
    y   = y(id).';
    gy  = gy(id,:);
    ggy = ggy(id,:,:);
    
    % Make a precision matrix
    W   = diag(w(id));

    % Build a matrix that Picks out the measured components of x
    I   = eye(nv);  % for chain rule
    PP  = I(id,:);  % picking/pick out the measured ones
    e   = PP*x - y; % calculated - measured (minus)
    ge  = PP - gy;
    
    nlam    = size(M,1) + size(K,1);
    lam     = z(nv+1:end);  % Lagrange multipliers 

    % constraint equations
    c   = [  M * q( x ); ...
             K * x ] ;
    f   = 0.5 * e.' * W * e  + lam.' * c ;  % limp, method of lagrange multipliers    
    % -(-1/2 sum of squares) + constraint eqns, minimize f => grad(f) = 0
    dcdx = [ M * diag( sys.dqdx( x ) ); ...
             K  ];
 
    g    = [ e.' * W * ge + lam.' * dcdx,  c.' ];

    ddq     =  diag( sys.d2qdx2( x ) ); % q"
        
    [nr,nc] = size(M);
    gg      = zeros(nc);
    for row = (1:nr)
        gg  = gg + lam(row)*diag(M(row,:))*ddq;
    end
    tmp = zeros(length(x));
    eW  = e.'*W;
    for jj = 1:size(ggy,1)
        tmp = tmp + eW(jj)*(squeeze(ggy(jj,:,:)));
    end
    H   = [  ge.'*W*ge-tmp+gg,  dcdx.'    ; ... % derivatives wrt lambdas
             dcdx            ,  zeros(nlam)  ]; % derivatives wrt var's
    g = g(:); % make sure g is returned as a column vector
    % keyboard
end

% -----------------------------------------------------------------------------------

function [opt] = check_opt(opt)
    % check opt input
    isbad = @(thing) (isempty(thing) & sum(isnan(thing)));

    % opt.printmes
    if ~isfield(opt,'printmes') || isbad(opt.printmes)
        opt.printmes = 1; % default on
    end
    % opt.K1K2
    if ~isfield(opt,'K1K2') || isbad(opt.K1K2)
        if opt.printmes ~= 0
            fprintf('No K1K2 formulation chosen. Assuming opt.K1K2 = 4\n');
        end
        opt.K1K2 = 4; % default K1K2 setting
    elseif opt.K1K2 > 18 || opt.K1K2 < 1 || ...
                opt.K1K2 == 6 || opt.K1K2 == 7 || opt.K1K2 == 8
        if opt.printmes ~= 0
            fprintf(['Invalid K1K2 formulation chosen. ' ...
                'Assuming opt.K1K2 = 4\n']);
        end
        opt.K1K2 = 4; % default K1K2 setting
    end
    % opt.TB
    if ~isfield(opt,'TB') || isbad(opt.TB)
        if opt.printmes ~= 0
            fprintf('No TB formulation chosen. Assuming opt.TB = 2\n');
        end
        opt.TB = 2;
    elseif opt.TB > 2 || opt.TB < 1
        if opt.printmes ~= 0
            fprintf(['Invalid TB formulation chosen. ' ...
                     'Assuming opt.TB = 2\n']);
        end
        opt.TB = 2;
    end
    % opt.KSO4
    if ~isfield(opt,'KSO4') || isbad(opt.KSO4)
        if opt.printmes ~= 0
            fprintf('No KSO4 formulation chosen. Assuming opt.KSO4 = 1\n');
        end
        opt.KSO4 = 1; % default opt.KSO4 setting
    elseif opt.KSO4 > 3 || opt.KSO4 < 1
        if opt.printmes ~= 0
            fprintf(['Invalid KSO4 formulation chosen. ' ...
                     'Assuming opt.KSO4 = 1\n']);
        end
        opt.KSO4 = 1; % default opt.KSO4 setting
    end
    % opt.KF
    if ~isfield(opt,'KF') || isbad(opt.KF)
        if opt.printmes ~= 0
            fprintf('No KF formulation chosen. Assuming opt.KF = 2 \n');
        end
        opt.KF = 2; % default KF
    elseif opt.KF > 2 || opt.KF < 1 
        if opt.printmes ~= 0
            fprintf(['Invalid KF formulation chosen. ' ...
                     'Assuming opt.KF = 2 \n']);
        end
        opt.KF = 2;
    end
    % opt.phscale
    if ~isfield(opt,'phscale') || isbad(opt.phscale)
        error(['No opt.phscale chosen, must choose 1 = tot, ' ...
               '2 = sws, 3 = free, 4 = NBS \n'])
    elseif opt.phscale > 4 || opt.phscale < 1
        eror(['Invalid opt.phscale chosen, must choose 1 = tot, ' ...
              '2 = sws, 3 = free, 4 = NBS \n'])
    end
    % opt.printcsv and opt.fname
    if ~isfield(opt,'printcsv') || isbad(opt.printcsv)
        opt.printcsv = 0; % default off
    elseif opt.printcsv > 1 || opt.printcsv < 0
        if opt.printmes ~= 0
            fprintf('Invalid CSV opt chosen. Assuming opt.csv = 1\n');
        end
    else
        if ~isfield(opt,'fname') || isbad(opt.fname)
            opt.fname = 'QUODcarb_output.csv';
            if opt.printmes ~= 0
                fprintf(['Invalid CSV filename. Assuming opt.fname' ...
                    ' = ''QUODcarb_output.csv'' \n']);
            end
        end
    end
    % opt.co2press
    if ~isfield(opt,'co2press') || isbad(opt.co2press)
        opt.co2press = 1; % on
        if opt.printmes ~=0
            fprintf('No opt.co2press chosen. Assuming opt.co2press = 1 (on). \n');
        end
    end
    % opt.Revelle
    if ~isfield(opt,'Revelle') || isbad(opt.Revelle)
        opt.Revelle = 0;
        if opt.printmes ~= 0
            fprintf('No opt.Revelle chosen. Assuming opt.Revelle = 0 (off). \n');
        end
    end
    % opt.turnoff totals (only got TB to work)
    if ~isfield(opt.turnoff,'TB') || isbad(opt.turnoff.TB)
        opt.turnoff.TB = 0; % default = not turned off
    end
    if ~isfield(opt.turnoff,'pK1') || isbad(opt.turnoff.pK1)
        opt.turnoff.pK1 = 0; % default = not turned off
    end
    % organic alkalinity
    if ~isfield(opt,'pKalpha') || isbad(opt.pKalpha)
        opt.pKalpha = 0; % off
    end
    if ~isfield(opt,'pKbeta') || isbad(opt.pKbeta)
        opt.pKbeta = 0; % off
    end
end

% ------------------------------------------------------------------------

function [y,gy,ggy] = update_y(y,x,obs,sys,opt)
    nTP         = length(sys.tp);
    M           = sys.M;
    K           = sys.K;
    nv          = size(M,2);
    gy          = zeros(nv,length(y));
    ggy         = zeros(nv,length(y),length(y));
    sal         = x(sys.isal);
    e3          = eps^3;
    ie3         = sqrt(-1)*e3;
    % calculate totals (see Ref's in calc_pTOT)
    [pT,gpT,ggpT,~] = calc_pTOT(opt,sal);
    pTB  = pT(1);  gpTB = gpT(1);  ggpTB = ggpT(1);
    pTS  = pT(2);  gpTS = gpT(2);  ggpTS = ggpT(2);
    pTF  = pT(3);  gpTF = gpT(3);  ggpTF = ggpT(3);
    pTCa = pT(4); gpTCa = gpT(4); ggpTCa = ggpT(4);
    % update y
    if (isnan(obs.TB))
        if opt.turnoff.TB ~= 1
            y(sys.ipTB)                     = pTB;
            gy(sys.ipTB,sys.isal)           = gpTB;
            ggy(sys.ipTB,sys.isal,sys.isal) = ggpTB;
        end
    end
    if (isnan(obs.TS))
        y(sys.ipTS)                     = pTS;
        gy(sys.ipTS,sys.isal)           = gpTS;
        ggy(sys.ipTS,sys.isal,sys.isal) = ggpTS;
    end
    if (isnan(obs.TF))
        y(sys.ipTF)                     = pTF;
        gy(sys.ipTF,sys.isal)           = gpTF;
        ggy(sys.ipTF,sys.isal,sys.isal) = ggpTF;
    end
    if (isnan(obs.TCa))
        y(sys.ipTCa)                        = pTCa;
        gy(sys.ipTCa,sys.isal)              = gpTCa;
        ggy(sys.ipTCa,sys.isal,sys.isal)    = ggpTCa;
    end

    for i = 1:nTP
        % use complex step method to get ∂T, ∂S, ∂P
        [pK,gpK] = calc_pK(opt,x(sys.tp(i).iT), x(sys.isal), x(sys.tp(i).iP));

        [pK_T,gpK_T] = calc_pK(opt, x(sys.tp(i).iT) + ie3  , x(sys.isal)       , x(sys.tp(i).iP)       );
        [pK_S,gpK_S] = calc_pK(opt, x(sys.tp(i).iT)        , x(sys.isal) + ie3 , x(sys.tp(i).iP)       );
        [pK_P,gpK_P] = calc_pK(opt, x(sys.tp(i).iT)        , x(sys.isal)       , x(sys.tp(i).iP) + ie3 );
        %
        ggpK = zeros(length(pK),3,3);
        ggpK(:,1,:) = imag(gpK_T)/e3;
        ggpK(:,2,:) = imag(gpK_S)/e3;
        ggpK(:,3,:) = imag(gpK_P)/e3;
        
        iTSP        = [ sys.tp(i).iT, sys.isal, sys.tp(i).iP];
        if (isnan(obs.tp(i).pK0))
            y(sys.tp(i).ipK0)             = pK(1); 
            gy(sys.tp(i).ipK0,iTSP)       = gpK(1,:);
            ggy(sys.tp(i).ipK0,iTSP,iTSP) = ggpK(1,:,:);
        end
        if (isnan(obs.tp(i).pK1))
            if opt.turnoff.pK1 ~= 1
                y(sys.tp(i).ipK1)             = pK(2);
                gy(sys.tp(i).ipK1,iTSP)       = gpK(2,:);
                ggy(sys.tp(i).ipK1,iTSP,iTSP) = ggpK(2,:,:);
            end
        end
        if (isnan(obs.tp(i).pK2))
            y(sys.tp(i).ipK2)             = pK(3);
            gy(sys.tp(i).ipK2,iTSP)       = gpK(3,:);
            ggy(sys.tp(i).ipK2,iTSP,iTSP) = ggpK(3,:,:);
        end   
        if (isnan(obs.tp(i).pKb))
            y(sys.tp(i).ipKb)             = pK(4);
            gy(sys.tp(i).ipKb,iTSP)       = gpK(4,:);
            ggy(sys.tp(i).ipKb,iTSP,iTSP) = ggpK(4,:,:);
        end
        if (isnan(obs.tp(i).pKw))
            y(sys.tp(i).ipKw)             = pK(5);
            gy(sys.tp(i).ipKw,iTSP)       = gpK(5,:);
            ggy(sys.tp(i).ipKw,iTSP,iTSP) = ggpK(5,:,:);
        end
        if (isnan(obs.tp(i).pKs))
            y(sys.tp(i).ipKs)             = pK(6);
            gy(sys.tp(i).ipKs,iTSP)       = gpK(6,:);
            ggy(sys.tp(i).ipKs,iTSP,iTSP) = ggpK(6,:,:);
        end
        if (isnan(obs.tp(i).pKf))
            y(sys.tp(i).ipKf)             = pK(7);
            gy(sys.tp(i).ipKf,iTSP)       = gpK(7,:);
            ggy(sys.tp(i).ipKf,iTSP,iTSP) = ggpK(7,:,:);
        end
        if (isnan(obs.tp(i).pKp1))
            y(sys.tp(i).ipKp1)             = pK(8); 
            gy(sys.tp(i).ipKp1,iTSP)       = gpK(8,:);
            ggy(sys.tp(i).ipKp1,iTSP,iTSP) = ggpK(8,:,:);
        end
        if (isnan(obs.tp(i).pKp2))
            y(sys.tp(i).ipKp2)             = pK(9);
            gy(sys.tp(i).ipKp2,iTSP)       = gpK(9,:);
            ggy(sys.tp(i).ipKp2,iTSP,iTSP) = ggpK(9,:,:);
        end
        if (isnan(obs.tp(i).pKp3))
            y(sys.tp(i).ipKp3)             = pK(10); 
            gy(sys.tp(i).ipKp3,iTSP)       = gpK(10,:);
            ggy(sys.tp(i).ipKp3,iTSP,iTSP) = ggpK(10,:,:);
        end  
        if (isnan(obs.tp(i).pKsi))
            y(sys.tp(i).ipKsi)             = pK(11);  
            gy(sys.tp(i).ipKsi,iTSP)       = gpK(11,:);
            ggy(sys.tp(i).ipKsi,iTSP,iTSP) = ggpK(11,:,:);
        end        
        if (isnan(obs.tp(i).pKnh4))
            y(sys.tp(i).ipKnh4)             = pK(12); 
            gy(sys.tp(i).ipKnh4,iTSP)       = gpK(12,:);
            ggy(sys.tp(i).ipKnh4,iTSP,iTSP) = ggpK(12,:,:);
        end   
        if (isnan(obs.tp(i).pKh2s))
            y(sys.tp(i).ipKh2s)             = pK(13); 
            gy(sys.tp(i).ipKh2s,iTSP)       = gpK(13,:); 
            ggy(sys.tp(i).ipKh2s,iTSP,iTSP) = ggpK(13,:,:); 
        end
        if (isnan(obs.tp(i).pp2f))
            y(sys.tp(i).ipp2f)             = pK(14);
            gy(sys.tp(i).ipp2f,iTSP)       = gpK(14,:);
            ggy(sys.tp(i).ipp2f,iTSP,iTSP) = ggpK(14,:,:);
        end
        if (isnan(obs.tp(i).pKar))
            y(sys.tp(i).ipKar)             = pK(15);
            gy(sys.tp(i).ipKar,iTSP)       = gpK(15,:);
            ggy(sys.tp(i).ipKar,iTSP,iTSP) = ggpK(15,:,:);
        end
        if (isnan(obs.tp(i).pKca))
            y(sys.tp(i).ipKca)             = pK(16);
            gy(sys.tp(i).ipKca,iTSP)       = gpK(16,:);
            ggy(sys.tp(i).ipKca,iTSP,iTSP) = ggpK(16,:,:);
        end
    end
end

% ------------------------------------------------------------------------

function z0 = init(opt,yobs,sys)
    q = sys.q;
    p = sys.p;
    
    y0  = yobs; 
    dic = q(yobs(sys.ipTC));
    alk = q(yobs(sys.ipTA));
    if (isnan(dic))
        dic         = 2200e-6;
        y0(sys.ipTC) = p(dic);
    end
    if (isnan(alk))
        alk         = 2200e-6;
        y0(sys.ipTA) = p(alk);
    end
    % calculate totals
    S = yobs(sys.isal);
    [pT,~,~,~] = calc_pTOT(opt,S);
    pTB  = pT(1);

    % calculate tp-independent vars
    gam     = dic/alk;

    if opt.turnoff.TB == 1
        TB = q(pTB);
        y0(sys.ipTB) = pTB;
    else
        TB = q(yobs(sys.ipTB));
        y0(sys.ipTB) = p(TB);
    end

    TS = q(yobs(sys.ipTS));
    y0(sys.ipTS) = p(TS);

    TF = q(yobs(sys.ipTF));
    y0(sys.ipTF) = p(TF);

    TP      = q(yobs(sys.ipTP));
    y0(sys.ipTP) = p(TP);

    TSi     = q(yobs(sys.ipTSi));
    y0(sys.ipTSi) = p(TSi);

    TNH4    = q(yobs(sys.ipTNH4));
    y0(sys.ipTNH4) = p(TNH4);

    TH2S    = q(yobs(sys.ipTH2S));
    y0(sys.ipTH2S) = p(TH2S);

    TCa     = q(yobs(sys.ipTCa));
    y0(sys.ipTCa) = p(TCa);

    if opt.pKalpha == 1 
        TAlpha = q(yobs(sys.ipTAlpha));
        y0(sys.ipTAlpha) = p(TAlpha);
    end
    if opt.pKbeta == 1
        TBeta = q(yobs(sys.ipTBeta));
        y0(sys.ipTBeta) = p(TBeta);
    end

    nTP = length(sys.tp);
    for i = 1:nTP
        % calculate pK1
        T = y0(sys.tp(i).iT); P = y0(sys.tp(i).iP);
        [pK,~,~] = calc_pK(opt,T,S,P);
        if opt.turnoff.pK1 == 1
            pK1 = pK(2);
            K1 = q(pK1);
            y0(sys.tp(i).ipK1) = pK1;
        else
            K1      = q(y0(sys.tp(i).ipK1));
        end

        % solve for the [H+] using only the carbonate alkalinity
        K0      = q(y0(sys.tp(i).ipK0));
        % K1      = q(y0(sys.tp(i).ipK1));
        K2      = q(y0(sys.tp(i).ipK2));
        h       = 0.5*( ( gam - 1 ) * K1 + ( ( 1 - gam )^2 * ...
                                             K1^2 - 4 * K1 * K2 * ( 1 - 2 * gam ) ).^0.5 ) ;
        hco3    = h * alk / (h + 2 * K2 );
        co2st   = h * hco3 / K1 ;
        co3     = dic*K1*K2/(K1*h + h*h + K1*K2) ;
        fco2    = co2st/K0;
        
        y0(sys.tp(i).iph)     = p(h);
        y0(sys.tp(i).iphco3)  = p(hco3);
        y0(sys.tp(i).ipco2st) = p(co2st);
        y0(sys.tp(i).ipco3)   = p(co3);
        y0(sys.tp(i).ipfco2)  = p(fco2);
        
        Kb      = q(y0(sys.tp(i).ipKb));
        boh4    = TB * Kb / (Kb + h) ;
        boh3    = TB - boh4;
        y0(sys.tp(i).ipboh3) = p(boh3);
        y0(sys.tp(i).ipboh4) = p(boh4);
        
        Kw      = q(y0(sys.tp(i).ipKw));
        oh      = Kw / h;
        y0(sys.tp(i).ipoh)   = p(oh);
        
        Ks      = q(y0(sys.tp(i).ipKs));
        h_tot   = h * ( 1 + TS / Ks );
        h_free  = h;
        hso4    = h_tot - h_free;
        so4     = Ks * hso4 / h;
        y0(sys.tp(i).iphso4)    = p(hso4);
        y0(sys.tp(i).ipso4)     = p(so4);
        
        Kf      = q(y0(sys.tp(i).ipKf));
        HF      = TF / ( 1 + Kf / h_free );
        F       = Kf * HF / h_free;
        h_sws   = h_tot + HF;
        y0(sys.tp(i).ipF)    = p(F);
        y0(sys.tp(i).ipHF)   = p(HF);
        
        Kp1     = q(y0(sys.tp(i).ipKp1));
        Kp2     = q(y0(sys.tp(i).ipKp2));
        Kp3     = q(y0(sys.tp(i).ipKp3));
        d       = ( h^3 + Kp1 * h^2 + Kp1 * Kp2 * h + Kp1 * Kp2 * Kp3);
        h3po4   = TP * h^3 / d;
        h2po4   = TP * Kp1 * h^2 / d;
        hpo4    = TP * Kp1 * Kp2 * h / d;
        po4     = TP * Kp1 * Kp2 * Kp3 / d;
        y0(sys.tp(i).iph3po4)    = p(h3po4);
        y0(sys.tp(i).iph2po4)    = p(h2po4);
        y0(sys.tp(i).iphpo4)     = p(hpo4);
        y0(sys.tp(i).ippo4)      = p(po4);
        
        Ksi     = q(y0(sys.tp(i).ipKsi));
        siooh3  = TSi / ( 1 + h / Ksi );
        sioh4   = TSi - siooh3;
        y0(sys.tp(i).ipsiooh3)   = p(siooh3);
        y0(sys.tp(i).ipsioh4)    = p(sioh4);
        
        Knh4    = q(y0(sys.tp(i).ipKnh4));
        nh3     = TNH4 / ( 1 + h / Knh4 );
        nh4     = TNH4 - nh3 ;
        y0(sys.tp(i).ipnh3)  = p(nh3);
        y0(sys.tp(i).ipnh4)  = p(nh4);
        
        Kh2s    = q(y0(sys.tp(i).ipKh2s));
        hs      = TH2S / ( 1 + h / Kh2s );
        h2s     = TH2S - hs ;
        y0(sys.tp(i).ipHS)   = p(hs);
        y0(sys.tp(i).ipH2S)  = p(h2s);
        
        p2f     = q(y0(sys.tp(i).ipp2f));
        pco2    = fco2/p2f;
        y0(sys.tp(i).ippco2)  = p(pco2);
        
        Kar     = q(y0(sys.tp(i).ipKar));
        OmegaAr = co3 * TCa / Kar;
        Kca     = q(y0(sys.tp(i).ipKca));
        OmegaCa = co3 * TCa / Kca ;
        y0(sys.tp(i).ipca)       = p(TCa);
        y0(sys.tp(i).ipOmegaAr)  = p(OmegaAr);
        y0(sys.tp(i).ipOmegaCa)  = p(OmegaCa);
                
        y0(sys.tp(i).iph_tot)   = p(h_tot);
        y0(sys.tp(i).iph_sws)   = p(h_sws);
        y0(sys.tp(i).iph_free)  = p(h_free);
        fH      = q(y0(sys.tp(i).ipfH));
        h_nbs   = h_free*fH;
        y0(sys.tp(i).iph_nbs)   = p(h_nbs);

        if opt.pKalpha == 1
            Kalpha = q(y0(sys.tp(i).ipKalpha));   
            alpha = (Kalpha*TAlpha)/(h+Kalpha);
            halpha = TAlpha - alpha;
            % alpha = (1/2)*Talpha; % make a function of pH
            % halpha = (1/2)*Talpha;
            y0(sys.tp(i).ipalpha) = p(alpha);
            y0(sys.tp(i).iphalpha) = p(halpha);
        end
        if opt.pKbeta == 1
            Kbeta = q(y0(sys.tp(i).ipKbeta));
            beta = (Kbeta*TBeta)/(h+Kbeta);
            hbeta = TBeta - beta;
            % beta = (1/2)*Tbeta;
            % hbeta = (1/2)*Tbeta;
            y0(sys.tp(i).ipbeta) = p(beta);
            y0(sys.tp(i).iphbeta) = p(hbeta);
        end
    end
    nlam    = size(sys.M,1) + size(sys.K,1);
    lam     = zeros(nlam,1);
    z0      = [y0(:);lam(:)];
end

% ------------------------------------------------------------------------

function [est] = parse_output(z,sigx,sys,f) % ORG ALK
    % populate est, output structure with best estimates
    %
    % INPUT:
    %   z  := system state including equilibrium constants and lagrange multipliers
    
    p       = sys.p;
    q       = sys.q;

    ebar    = @(j) (0.5 * ( q( z(j) - sigx(j) ) - q( z(j) + sigx(j) ) ) );
    ebar_l  = @(j) ( q( -sigx(j) ) ); % lower sigma
    ebar_u  = @(j) ( q( sigx(j) ) ); % upper sigma
        
    % populate 'est' structure with best estimate:
    %   1. p(value) and p(error) where p(x) = -log10(x)
    %   2. value and average error about the value in 'q' 
    %           where q(x) = x^(-10)
    %   3. upper and lower bounds in 'q' space, not symmetric
    %           about the value in 'q' space
    
    est.f       = f; % residual f value, from limp

    est.sal     = z(sys.isal);
    est.esal    = sigx(sys.isal);
    
    % TC (DIC)
    est.pTC     = z(sys.ipTC);               
    est.epTC    = sigx(sys.ipTC);    
    est.TC      = q(z(sys.ipTC))*1e6; % 1e6  converts mol/kg to µmol/kg
    est.eTC     = ebar(sys.ipTC)*1e6;      
    est.eTC_l   = ebar_l(sys.ipTC)*1e6;    
    est.eTC_u   = ebar_u(sys.ipTC)*1e6;
    
    % TA Alkalinity
    est.pTA     = z(sys.ipTA);               
    est.epTA    = sigx(sys.ipTA);
    est.TA      = q(z(sys.ipTA))*1e6; % convt 
    est.eTA     = ebar(sys.ipTA)*1e6;      
    est.eTA_l   = ebar_l(sys.ipTA)*1e6;    
    est.eTA_u   = ebar_u(sys.ipTA)*1e6;

    % TB borate
    est.pTB     = z(sys.ipTB);
    est.epTB    = sigx(sys.ipTB);
    est.TB      = q(z(sys.ipTB))*1e6; % convt mol/kg to µmol/kg
    est.eTB     = ebar(sys.ipTB)*1e6;
    est.eTB_l   = ebar_l(sys.ipTB)*1e6;
    est.eTB_u   = ebar_u(sys.ipTB)*1e6;
    
    % TS sulfate
    est.pTS     = z(sys.ipTS);
    est.epTS    = sigx(sys.ipTS);
    est.TS      = q(z(sys.ipTS))*1e6;  % convt mol/kg to µmol/kg
    est.eTS     = ebar(sys.ipTS)*1e6;
    est.eTS_l   = ebar_l(sys.ipTS)*1e6;
    est.eTS_u   = ebar_u(sys.ipTS)*1e6;
    
    % TF fluoride
    est.pTF     = z(sys.ipTF);
    est.epTF        = sigx(sys.ipTF);
    est.TF      = q(z(sys.ipTF))*1e6; % convt mol/kg to µmol/kg
    est.eTF     = ebar(sys.ipTF)*1e6;
    est.eTF_l   = ebar_l(sys.ipTF)*1e6;
    est.eTF_u   = ebar_u(sys.ipTF)*1e6;

    % TP Phosphate
    est.pTP     = z(sys.ipTP);           
    est.epTP    = sigx(sys.ipTP);
    est.TP      = q(z(sys.ipTP))*1e6;   % convt mol/kg to µmol/kg   
    est.eTP     = ebar(sys.ipTP)*1e6;
    est.eTP_l   = ebar_l(sys.ipTP)*1e6; 
    est.eTP_u   = ebar_u(sys.ipTP)*1e6;

    % TSi silicate
    est.pTSi    = z(sys.ipTSi);         
    est.epTSi   = sigx(sys.ipTSi);
    est.TSi     = q(z(sys.ipTSi))*1e6;  % convt mol/kg to µmol/kg 
    est.eTSi    = ebar(sys.ipTSi)*1e6;
    est.eTSi_l  = ebar_l(sys.ipTSi)*1e6; 
    est.eTSi_u  = ebar_u(sys.ipTSi)*1e6;
    
    % TNH4 nitrate
    est.pTNH4       = z(sys.ipTNH4);       
    est.epTNH4      = sigx(sys.ipTNH4);
    est.TNH4        = q(z(sys.ipTNH4))*1e6;   % convt mol/kg to µmol/kg
    est.eTNH4       = ebar(sys.ipTNH4)*1e6;
    est.eTNH4_l     = ebar_l(sys.ipTNH4)*1e6; 
    est.eTNH4_u     = ebar_u(sys.ipTNH4)*1e6;
    
    % TH2S sulfide
    est.pTH2S       = z(sys.ipTH2S);       
    est.epTH2S      = sigx(sys.ipTH2S);
    est.TH2S        = q(z(sys.ipTH2S))*1e6; % convt mol/kg to µmol/kg
    est.eTH2S       = ebar(sys.ipTH2S)*1e6;
    est.eTH2S_l     = ebar_l(sys.ipTH2S)*1e6; 
    est.eTH2S_u     = ebar_u(sys.ipTH2S)*1e6;

    % TCa calcium
    est.pTCa       = z(sys.ipTCa);       
    est.epTCa      = sigx(sys.ipTCa);
    est.TCa        = q(z(sys.ipTCa))*1e6; % convt mol/kg to µmol/kg
    est.eTCa       = ebar(sys.ipTCa)*1e6;
    est.eTCa_l     = ebar_l(sys.ipTCa)*1e6; 
    est.eTCa_u     = ebar_u(sys.ipTCa)*1e6;
    
    if (isfield(sys,'ipTAlpha'))
        est.pTAlpha     = z(sys.ipTAlpha);
        est.epTAlpha    = sigx(sys.ipTAlpha);
        est.TAlpha      = q(z(sys.ipTAlpha))*1e6;
        est.eTAlpha     = ebar(sys.ipTAlpha)*1e6;
    end
    if (isfield(sys,'ipTBeta'))
        est.pTBeta      = z(sys.ipTBeta);
        est.epTBeta     = sigx(sys.ipTBeta);
        est.TBeta       = q(z(sys.ipTBeta))*1e6;
        est.eTBeta      = ebar(sys.ipTBeta)*1e6;
    end

    nTP = length(sys.tp);
    for i = 1:nTP
        % temp (deg C)
        est.tp(i).T     = z(sys.tp(i).iT);
        est.tp(i).eT    = sigx(sys.tp(i).iT);
        est.tp(i).eT_l  = z(sys.tp(i).iT)-sigx(sys.tp(i).iT);
        est.tp(i).eT_u  = z(sys.tp(i).iT)+sigx(sys.tp(i).iT);
        
        % pressure (dbar)
        est.tp(i).P     = z(sys.tp(i).iP);        
        est.tp(i).eP    = sigx(sys.tp(i).iP);
        est.tp(i).eP_l  = z(sys.tp(i).iP)-sigx(sys.tp(i).iP);
        est.tp(i).eP_u  = z(sys.tp(i).iP)+sigx(sys.tp(i).iP);
        
        % fCO2
        est.tp(i).fco2      = q(z(sys.tp(i).ipfco2)) * 1e6; % convt atm to µatm
        est.tp(i).efco2     = ebar(sys.tp(i).ipfco2) * 1e6;
        est.tp(i).efco2_l   = ebar_l(sys.tp(i).ipfco2) * 1e6;
        est.tp(i).efco2_u   = ebar_u(sys.tp(i).ipfco2) * 1e6;
        est.tp(i).pfco2     = z(sys.tp(i).ipfco2);
        est.tp(i).epfco2    = sigx(sys.tp(i).ipfco2);

        % pCO2
        est.tp(i).pco2      = q(z(sys.tp(i).ippco2)) * 1e6; % convt atm to µatm
        est.tp(i).epco2     = ebar(sys.tp(i).ippco2) * 1e6;
        est.tp(i).epco2_l   = ebar_l(sys.tp(i).ippco2) * 1e6;
        est.tp(i).epco2_u   = ebar_u(sys.tp(i).ippco2) * 1e6;
        est.tp(i).ppco2     = z(sys.tp(i).ippco2);
        est.tp(i).eppco2    = sigx(sys.tp(i).ippco2);
        
        % HCO3
        est.tp(i).hco3      = q(z(sys.tp(i).iphco3))*1e6; % convt mol/kg to µmol/kg
        est.tp(i).ehco3     = ebar(sys.tp(i).iphco3)*1e6;
        est.tp(i).ehco3_l   = ebar_l(sys.tp(i).iphco3)*1e6;
        est.tp(i).ehco3_u   = ebar_u(sys.tp(i).iphco3)*1e6;
        est.tp(i).phco3     = z(sys.tp(i).iphco3);
        est.tp(i).ephco3    = sigx(sys.tp(i).iphco3);

        % CO3
        est.tp(i).co3       = q(z(sys.tp(i).ipco3))*1e6; % convt mol/kg to µmol/kg
        est.tp(i).eco3      = ebar(sys.tp(i).ipco3)*1e6;
        est.tp(i).eco3_l    = ebar_l(sys.tp(i).ipco3)*1e6;
        est.tp(i).eco3_u    = ebar_u(sys.tp(i).ipco3)*1e6;
        est.tp(i).pco3      = z(sys.tp(i).ipco3);
        est.tp(i).epco3     = sigx(sys.tp(i).ipco3);

        % CO2*
        est.tp(i).co2st     = q(z(sys.tp(i).ipco2st)); % convt mol/kg to µmol/kg
        est.tp(i).eco2st    = ebar(sys.tp(i).ipco2st);
        est.tp(i).eco2st_l  = ebar_l(sys.tp(i).ipco2st);
        est.tp(i).eco2st_u  = ebar_u(sys.tp(i).ipco2st);
        est.tp(i).pco2st    = z(sys.tp(i).ipco2st);
        est.tp(i).epco2st   = sigx(sys.tp(i).ipco2st);

        % pH on the scale opt.phscale used to compute the pK values
        est.tp(i).ph        = z(sys.tp(i).iph); % (9)
        est.tp(i).eph       = sigx(sys.tp(i).iph);
        est.tp(i).h         = q(z(sys.tp(i).iph)) * 1e6;
        est.tp(i).eh        = ebar(sys.tp(i).iph) * 1e6;
        est.tp(i).eh_l      = ebar_l(sys.tp(i).iph) * 1e6;
        est.tp(i).eh_u      = ebar_u(sys.tp(i).iph) * 1e6;

        % pH_free
        est.tp(i).ph_free    = z(sys.tp(i).iph_free); % (15)
        est.tp(i).eph_free   = sigx(sys.tp(i).iph_free);
        est.tp(i).h_free     = q(z(sys.tp(i).iph_free)) * 1e6;
        est.tp(i).eh_free    = ebar(sys.tp(i).iph_free) * 1e6;
        est.tp(i).eh_free_l  = ebar_l(sys.tp(i).iph_free) * 1e6;
        est.tp(i).eh_free_u  = ebar_u(sys.tp(i).iph_free) * 1e6;

        % pH_tot
        est.tp(i).ph_tot    = z(sys.tp(i).iph_tot); % (21)
        est.tp(i).eph_tot   = sigx(sys.tp(i).iph_tot);
        est.tp(i).h_tot     = q(z(sys.tp(i).iph_tot)) * 1e6;
        est.tp(i).eh_tot    = ebar(sys.tp(i).iph_tot) * 1e6;
        est.tp(i).eh_tot_l  = ebar_l(sys.tp(i).iph_tot) * 1e6;
        est.tp(i).eh_tot_u  = ebar_u(sys.tp(i).iph_tot) * 1e6;

        % pH_sws
        est.tp(i).ph_sws    = z(sys.tp(i).iph_sws);
        est.tp(i).eph_sws   = sigx(sys.tp(i).iph_sws);
        est.tp(i).h_sws     = q(z(sys.tp(i).iph_sws)) * 1e6;
        est.tp(i).eh_sws    = ebar(sys.tp(i).iph_sws) * 1e6;
        est.tp(i).eh_sws_l  = ebar_l(sys.tp(i).iph_sws) * 1e6;
        est.tp(i).eh_sws_u  = ebar_u(sys.tp(i).iph_sws) * 1e6;

        % pH_nbs
        est.tp(i).ph_nbs    = z(sys.tp(i).iph_nbs);
        est.tp(i).eph_nbs   = sigx(sys.tp(i).iph_nbs);
        est.tp(i).h_nbs     = q(z(sys.tp(i).iph_nbs)) * 1e6;
        est.tp(i).eh_nbs    = ebar(sys.tp(i).iph_nbs) * 1e6;
        est.tp(i).eh_nbs_l  = ebar_l(sys.tp(i).iph_nbs) * 1e6;
        est.tp(i).eh_nbs_u  = ebar_u(sys.tp(i).iph_nbs) * 1e6;

        % fH = activity coefficient
        est.tp(i).fH        = q(z(sys.tp(i).ipfH)) * 1e6;
        est.tp(i).efH       = ebar(sys.tp(i).ipfH) * 1e6;
        est.tp(i).efH_l     = ebar_l(sys.tp(i).ipfH) * 1e6;
        est.tp(i).efH_u     = ebar_u(sys.tp(i).ipfH) * 1e6;
        est.tp(i).pfH       = z(sys.tp(i).ipfH);
        est.tp(i).epfH      = sigx(sys.tp(i).ipfH);

        % p2f
        est.tp(i).p2f       = q(z(sys.tp(i).ipp2f));
        est.tp(i).ep2f      = ebar(sys.tp(i).ipp2f);
        est.tp(i).pp2f      = z(sys.tp(i).ipp2f); 
        est.tp(i).epp2f     = sigx(sys.tp(i).ipp2f);

        % pK0 
        est.tp(i).pK0       = z(sys.tp(i).ipK0); % (79)
        est.tp(i).epK0      = sigx(sys.tp(i).ipK0);
        est.tp(i).K0        = q(z(sys.tp(i).ipK0));
        est.tp(i).eK0       = ebar(sys.tp(i).ipK0);
        est.tp(i).eK0_l     = ebar_l(sys.tp(i).ipK0);
        est.tp(i).eK0_u     = ebar_u(sys.tp(i).ipK0);

        % pK1
        est.tp(i).pK1   = z(sys.tp(i).ipK1);
        est.tp(i).epK1  = sigx(sys.tp(i).ipK1);
        est.tp(i).K1    = q(z(sys.tp(i).ipK1));
        est.tp(i).eK1   = ebar(sys.tp(i).ipK1);
        est.tp(i).eK1_l = ebar_l(sys.tp(i).ipK1); 
        est.tp(i).eK1_u = ebar_u(sys.tp(i).ipK1);

        % pK2
        est.tp(i).pK2   = z(sys.tp(i).ipK2);
        est.tp(i).epK2  = sigx(sys.tp(i).ipK2);
        est.tp(i).K2    = q(z(sys.tp(i).ipK2));
        est.tp(i).eK2   = ebar(sys.tp(i).ipK2);
        est.tp(i).eK2_l = ebar_l(sys.tp(i).ipK2);
        est.tp(i).eK2_u = ebar_u(sys.tp(i).ipK2);
        
        % OH 
        est.tp(i).oh    = q(z(sys.tp(i).ipoh))*1e6; % convt
        est.tp(i).eoh   = ebar(sys.tp(i).ipoh)*1e6;
        est.tp(i).eoh_l = ebar_l(sys.tp(i).ipoh)*1e6;
        est.tp(i).eoh_u = ebar_u(sys.tp(i).ipoh)*1e6;
        est.tp(i).poh   = z(sys.tp(i).ipoh);
        est.tp(i).epoh  = sigx(sys.tp(i).ipoh);

        % pKw 
        est.tp(i).pKw   = z(sys.tp(i).ipKw); % (103)
        est.tp(i).epKw  = sigx(sys.tp(i).ipKw);
        est.tp(i).Kw    = q(z(sys.tp(i).ipKw));
        est.tp(i).eKw   = ebar(sys.tp(i).ipKw);
        est.tp(i).eKw_l = ebar_l(sys.tp(i).ipKw);
        est.tp(i).eKw_u = ebar_u(sys.tp(i).ipKw);

        % BOH4 borate
        est.tp(i).boh4      = q(z(sys.tp(i).ipboh4))*1e6; % convt mol/kg to µmol/kg
        est.tp(i).eboh4     = ebar(sys.tp(i).ipboh4)*1e6;
        est.tp(i).eboh4_l   = ebar_l(sys.tp(i).ipboh4)*1e6;
        est.tp(i).eboh4_u   = ebar_u(sys.tp(i).ipboh4)*1e6;
        est.tp(i).pboh4     = z(sys.tp(i).ipboh4);
        est.tp(i).epboh4    = sigx(sys.tp(i).ipboh4);

        % BOH3
        est.tp(i).boh3      = q(z(sys.tp(i).ipboh3))*1e6;
        est.tp(i).eboh3     = ebar(sys.tp(i).ipboh3)*1e6;
        est.tp(i).eboh3_l   = ebar_l(sys.tp(i).ipboh3)*1e6;
        est.tp(i).eboh3_u   = ebar_u(sys.tp(i).ipboh3)*1e6;
        est.tp(i).pboh3     = z(sys.tp(i).ipboh3);
        est.tp(i).epboh3    = sigx(sys.tp(i).ipboh3);
            
        % pKb
        est.tp(i).pKb       = z(sys.tp(i).ipKb); % (121)
        est.tp(i).epKb      = sigx(sys.tp(i).ipKb);
        est.tp(i).Kb        = q(z(sys.tp(i).ipKb));
        est.tp(i).eKb       = ebar(sys.tp(i).ipKb);
        est.tp(i).eKb_l     = ebar_l(sys.tp(i).ipKb);
        est.tp(i).eKb_u     = ebar_u(sys.tp(i).ipKb);

        % SO4 sulfate
        est.tp(i).so4       = q(z(sys.tp(i).ipso4)); % mol/kg
        est.tp(i).eso4      = ebar(sys.tp(i).ipso4);
        est.tp(i).eso4_l    = ebar_l(sys.tp(i).ipso4);
        est.tp(i).eso4_u    = ebar_u(sys.tp(i).ipso4);
        est.tp(i).pso4      = z(sys.tp(i).ipso4);
        est.tp(i).epso4     = sigx(sys.tp(i).ipso4);

        % HSO4
        est.tp(i).hso4      = q(z(sys.tp(i).iphso4))*1e6;
        est.tp(i).ehso4     = ebar(sys.tp(i).iphso4)*1e6;
        est.tp(i).ehso4_l   = ebar_l(sys.tp(i).iphso4)*1e6;
        est.tp(i).ehso4_u   = ebar_u(sys.tp(i).iphso4)*1e6;
        est.tp(i).phso4     = z(sys.tp(i).iphso4);
        est.tp(i).ephso4    = sigx(sys.tp(i).iphso4);

        % pKs
        est.tp(i).pKs       = z(sys.tp(i).ipKs); % (145)
        est.tp(i).epKs      = sigx(sys.tp(i).ipKs);
        est.tp(i).Ks        = q(z(sys.tp(i).ipKs));
        est.tp(i).eKs       = ebar(sys.tp(i).ipKs);
        est.tp(i).eKs_l     = ebar_l(sys.tp(i).ipKs);
        est.tp(i).eKs_u     = ebar_u(sys.tp(i).ipKs);

        % F fluoride 
        est.tp(i).F         = q(z(sys.tp(i).ipF))*1e6; % convt
        est.tp(i).eF        = ebar(sys.tp(i).ipF)*1e6;
        est.tp(i).eF_l      = ebar_l(sys.tp(i).ipF)*1e6;
        est.tp(i).ef_u      = ebar_u(sys.tp(i).ipF)*1e6;
        est.tp(i).pF        = z(sys.tp(i).ipF);
        est.tp(i).epF       = sigx(sys.tp(i).ipF);

        % HF 
        est.tp(i).HF        = q(z(sys.tp(i).ipHF))*1e6;
        est.tp(i).eHF       = ebar(sys.tp(i).ipHF)*1e6;
        est.tp(i).eHF_l     = ebar_l(sys.tp(i).ipHF)*1e6;
        est.tp(i).eHF_u     = ebar_u(sys.tp(i).ipHF)*1e6;
        est.tp(i).pHF       = z(sys.tp(i).ipHF);
        est.tp(i).epHF      = sigx(sys.tp(i).ipHF);

        % pKf
        est.tp(i).pKf       = z(sys.tp(i).ipKf); % (163)
        est.tp(i).epKf      = sigx(sys.tp(i).ipKf);
        est.tp(i).Kf        = q(z(sys.tp(i).ipKf));
        est.tp(i).eKf       = ebar(sys.tp(i).ipKf);
        est.tp(i).eKf_l     = ebar_l(sys.tp(i).ipKf);
        est.tp(i).eKf_u     = ebar_u(sys.tp(i).ipKf);

        % PO4
        est.tp(i).po4       = q(z(sys.tp(i).ippo4))*1e6; % convt
        est.tp(i).epo4      = ebar(sys.tp(i).ippo4)*1e6;
        est.tp(i).epo4_l    = ebar_l(sys.tp(i).ippo4)*1e6;
        est.tp(i).epo4_u    = ebar_u(sys.tp(i).ippo4)*1e6;
        est.tp(i).ppo4      = z(sys.tp(i).ippo4);
        est.tp(i).eppo4     = sigx(sys.tp(i).ippo4);

        % HPO4
        est.tp(i).hpo4      = q(z(sys.tp(i).iphpo4))*1e6;
        est.tp(i).ehpo4     = ebar(sys.tp(i).iphpo4)*1e6;
        est.tp(i).ehpo4_l   = ebar_l(sys.tp(i).iphpo4)*1e6;
        est.tp(i).ehpo4_u   = ebar_u(sys.tp(i).iphpo4)*1e6;
        est.tp(i).phpo4     = z(sys.tp(i).iphpo4);
        est.tp(i).ephpo4    = sigx(sys.tp(i).iphpo4);

        % H2PO4
        est.tp(i).h2po4     = q(z(sys.tp(i).iph2po4))*1e6;
        est.tp(i).eh2po4    = ebar(sys.tp(i).iph2po4)*1e6;
        est.tp(i).eh2po4_l  = ebar_l(sys.tp(i).iph2po4)*1e6;
        est.tp(i).eh2po4_u  = ebar_u(sys.tp(i).iph2po4)*1e6;
        est.tp(i).ph2po4    = z(sys.tp(i).iph2po4);
        est.tp(i).eph2po4   = sigx(sys.tp(i).iph2po4);            

        % H3PO4
        est.tp(i).h3po4     = q(z(sys.tp(i).iph3po4))*1e6;
        est.tp(i).eh3po4    = ebar(sys.tp(i).iph3po4)*1e6;
        est.tp(i).eh3po4_l  = ebar_l(sys.tp(i).iph3po4)*1e6;
        est.tp(i).eh3po4_u  = ebar_u(sys.tp(i).iph3po4)*1e6;
        est.tp(i).ph3po4    = z(sys.tp(i).iph3po4);
        est.tp(i).eph3po4   = sigx(sys.tp(i).iph3po4);

        % pKp1
        est.tp(i).pKp1      = z(sys.tp(i).ipKp1);
        est.tp(i).epKp1     = sigx(sys.tp(i).ipKp1);
        est.tp(i).Kp1       = q(z(sys.tp(i).ipKp1));
        est.tp(i).eKp1      = ebar(sys.tp(i).ipKp1);
        est.tp(i).eKp1_l    = ebar_l(sys.tp(i).ipKp1);
        est.tp(i).eKp1_u    = ebar_u(sys.tp(i).ipKp1);

        % pKp2
        est.tp(i).pKp2      = z(sys.tp(i).ipKp2);
        est.tp(i).epKp2     = sigx(sys.tp(i).ipKp2);
        est.tp(i).Kp2       = q(z(sys.tp(i).ipKp2));
        est.tp(i).eKp2      = ebar(sys.tp(i).ipKp2);
        est.tp(i).pKp2_l    = ebar_l(sys.tp(i).ipKp2);
        est.tp(i).eKp2_u    = ebar_u(sys.tp(i).ipKp2);

        % pKp3
        est.tp(i).pKp3      = z(sys.tp(i).ipKp3);
        est.tp(i).epKp3     = sigx(sys.tp(i).ipKp3);
        est.tp(i).Kp3       = q(z(sys.tp(i).ipKp3));
        est.tp(i).eKp3      = ebar(sys.tp(i).ipKp3);
        est.tp(i).eKp3_l    = ebar_l(sys.tp(i).ipKp3);
        est.tp(i).eKp3_u    = ebar_u(sys.tp(i).ipKp3);
        
        % SiOH4
        est.tp(i).sioh4     = q(z(sys.tp(i).ipsioh4))*1e6; % convt
        est.tp(i).esioh4    = ebar(sys.tp(i).ipsioh4)*1e6;
        est.tp(i).esioh4_l  = ebar_l(sys.tp(i).ipsioh4)*1e6;
        est.tp(i).esioh4_u  = ebar_u(sys.tp(i).ipsioh4)*1e6;
        est.tp(i).psioh4    = z(sys.tp(i).ipsioh4);
        est.tp(i).epsioh4   = sigx(sys.tp(i).ipsioh4);

        % SiOH3
        est.tp(i).siooh3    = q(z(sys.tp(i).ipsiooh3))*1e6;
        est.tp(i).esiooh3   = ebar(sys.tp(i).ipsiooh3)*1e6;
        est.tp(i).esiooh3_l = ebar_l(sys.tp(i).ipsiooh3)*1e6;
        est.tp(i).esiooh3_u = ebar_u(sys.tp(i).ipsiooh3)*1e6;
        est.tp(i).psiooh3   = z(sys.tp(i).ipsiooh3);
        est.tp(i).epsiooh3  = sigx(sys.tp(i).ipsiooh3);

        % pKsi
        est.tp(i).pKsi      = z(sys.tp(i).ipKsi);
        est.tp(i).epKsi     = sigx(sys.tp(i).ipKsi);
        est.tp(i).Ksi       = q(z(sys.tp(i).ipKsi));
        est.tp(i).eKsi      = ebar(sys.tp(i).ipKsi);
        est.tp(i).eKsi_l    = ebar_l(sys.tp(i).ipKsi);
        est.tp(i).eKsi_u    = ebar_u(sys.tp(i).ipKsi);

        % NH3
        est.tp(i).nh3       = q(z(sys.tp(i).ipnh3))*1e6; % convt
        est.tp(i).enh3      = ebar(sys.tp(i).ipnh3)*1e6;
        est.tp(i).enh3_l    = ebar_l(sys.tp(i).ipnh3)*1e6;
        est.tp(i).enh3_u    = ebar_u(sys.tp(i).ipnh3)*1e6;
        est.tp(i).pnh3      = z(sys.tp(i).ipnh3);
        est.tp(i).epnh3     = sigx(sys.tp(i).ipnh3);

        % NH4
        est.tp(i).nh4       = q(z(sys.tp(i).ipnh4))*1e6;
        est.tp(i).enh4      = ebar(sys.tp(i).ipnh4)*1e6;
        est.tp(i).enh4_l    = ebar_l(sys.tp(i).ipnh4)*1e6;
        est.tp(i).enh4_u    = ebar_u(sys.tp(i).ipnh4)*1e6;
        est.tp(i).pnh4      = z(sys.tp(i).ipnh4);
        est.tp(i).epnh4     = sigx(sys.tp(i).ipnh4);

        % pKNH4
        est.tp(i).pKnh4     = z(sys.tp(i).ipKnh4);
        est.tp(i).epKnh4    = sigx(sys.tp(i).ipKnh4);
        est.tp(i).Knh4      = q(z(sys.tp(i).ipKnh4));
        est.tp(i).eKnh4     = ebar(sys.tp(i).ipKnh4);
        est.tp(i).eKnh4_l   = ebar_l(sys.tp(i).ipKnh4);
        est.tp(i).eKnh4_u   = ebar_u(sys.tp(i).ipKnh4);

        % HS
        est.tp(i).HS        = q(z(sys.tp(i).ipHS))*1e6; % convt
        est.tp(i).eHS       = ebar(sys.tp(i).ipHS)*1e6;
        est.tp(i).eHS_l     = ebar_l(sys.tp(i).ipHS)*1e6;
        est.tp(i).eHS_u     = ebar_u(sys.tp(i).ipHS)*1e6;
        est.tp(i).pHS       = z(sys.tp(i).ipHS);
        est.tp(i).epHS      = sigx(sys.tp(i).ipHS);

        % H2S
        est.tp(i).H2S       = q(z(sys.tp(i).ipH2S))*1e6;
        est.tp(i).eH2S      = ebar(sys.tp(i).ipH2S)*1e6;
        est.tp(i).eH2S_l    = ebar_l(sys.tp(i).ipH2S)*1e6;
        est.tp(i).eHS2_u    = ebar_u(sys.tp(i).ipH2S)*1e6;
        est.tp(i).pH2S      = z(sys.tp(i).ipH2S);
        est.tp(i).epH2S     = sigx(sys.tp(i).ipH2S);

        % pKh2s
        est.tp(i).pKh2s     = z(sys.tp(i).ipKh2s);
        est.tp(i).epKh2s    = sigx(sys.tp(i).ipKh2s);
        est.tp(i).Kh2s      = q(z(sys.tp(i).ipKh2s));
        est.tp(i).eKh2s     = ebar(sys.tp(i).ipKh2s);
        est.tp(i).eKh2s_l   = ebar_l(sys.tp(i).ipKh2s);
        est.tp(i).eKh2s_u   = ebar_u(sys.tp(i).ipKh2s);

        % Ca
        est.tp(i).ca        = q(z(sys.tp(i).ipca))*1e6;
        est.tp(i).eca       = ebar(sys.tp(i).ipca)*1e6;
        est.tp(i).eca_l     = ebar_l(sys.tp(i).ipca)*1e6;
        est.tp(i).eca_u     = ebar_u(sys.tp(i).ipca)*1e6;
        est.tp(i).pca       = z(sys.tp(i).ipca);
        est.tp(i).epca      = sigx(sys.tp(i).ipca);

        % Omega_Ar
        est.tp(i).OmegaAr    = q(z(sys.tp(i).ipOmegaAr)); % unitless
        est.tp(i).eOmegaAr   = ebar(sys.tp(i).ipOmegaAr);
        est.tp(i).eOmegaAr_l = ebar_l(sys.tp(i).ipOmegaAr);
        est.tp(i).eOmegaAr_u = ebar_u(sys.tp(i).ipOmegaAr);
        est.tp(i).pOmegaAr   = z(sys.tp(i).ipOmegaAr);
        est.tp(i).epOmegaAr  = sigx(sys.tp(i).ipOmegaAr);

        % pKar
        est.tp(i).pKar      = z(sys.tp(i).ipKar);
        est.tp(i).epKar     = sigx(sys.tp(i).ipKar);
        est.tp(i).Kar       = q(z(sys.tp(i).ipKar));
        est.tp(i).eKar      = ebar(sys.tp(i).ipKar);
        est.tp(i).eKar_l    = ebar_l(sys.tp(i).ipKar);
        est.tp(i).eKar_u    = ebar_u(sys.tp(i).ipKar);

        % Omega_Ca
        est.tp(i).OmegaCa    = q(z(sys.tp(i).ipOmegaCa));
        est.tp(i).eOmegaCa   = ebar(sys.tp(i).ipOmegaCa);
        est.tp(i).eOmegaCa_l = ebar_l(sys.tp(i).ipOmegaCa);
        est.tp(i).eOmegaCa_u = ebar_u(sys.tp(i).ipOmegaCa);
        est.tp(i).pOmegaCa   = z(sys.tp(i).ipOmegaCa);
        est.tp(i).epOmegaCa  = sigx(sys.tp(i).ipOmegaCa);

        % pKca
        est.tp(i).pKca      = z(sys.tp(i).ipKca);
        est.tp(i).epKca     = sigx(sys.tp(i).ipKca);
        est.tp(i).Kca       = q(z(sys.tp(i).ipKca));
        est.tp(i).eKca      = ebar(sys.tp(i).ipKca);
        est.tp(i).eKca_l    = ebar_l(sys.tp(i).ipKca);
        est.tp(i).eKca_u    = ebar_u(sys.tp(i).ipKca);

        if (isfield(sys,'ipTAlpha'))
            % pKalpha
            est.tp(i).pKalpha   = z(sys.tp(i).ipKalpha);
            est.tp(i).epKalpha  = sigx(sys.tp(i).ipKalpha);
            % alpha
            est.tp(i).palpha    = z(sys.tp(i).ipalpha);
            est.tp(i).epalpha   = sigx(sys.tp(i).ipalpha);
            est.tp(i).alpha     = q(z(sys.tp(i).ipalpha))*1e6;
            est.tp(i).ealpha    = ebar(sys.tp(i).ipalpha)*1e6;
            % halpha
            est.tp(i).phalpha   = z(sys.tp(i).iphalpha);
            est.tp(i).ephalpha  = sigx(sys.tp(i).iphalpha);
            est.tp(i).halpha    = q(z(sys.tp(i).iphalpha))*1e6;
            est.tp(i).ehalpha   = ebar(sys.tp(i).iphalpha)*1e6;
        end

        if (isfield(sys,'ipTBeta'))
            % pKbeta
            est.tp(i).pKbeta = z(sys.tp(i).ipKbeta);
            est.tp(i).epKbeta = sigx(sys.tp(i).ipKbeta);
            % beta
            est.tp(i).pbeta = z(sys.tp(i).ipbeta);
            est.tp(i).epbeta = sigx(sys.tp(i).ipbeta);
            est.tp(i).beta = q(z(sys.tp(i).ipbeta))*1e6;
            est.tp(i).ebeta = ebar(sys.tp(i).ipbeta)*1e6;
            % hbeta
            est.tp(i).phbeta = z(sys.tp(i).iphbeta);
            est.tp(i).ephbeta = sigx(sys.tp(i).iphbeta);
            est.tp(i).hbeta = q(z(sys.tp(i).iphbeta))*1e6;
            est.tp(i).ehbeta = ebar(sys.tp(i).iphbeta)*1e6;
        end

    end
end





