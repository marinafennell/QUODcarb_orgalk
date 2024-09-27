function [obs,yobs,wobs,sys] = parse_input(obs,sys,opt,nD)
    % ORG ALK
    isgood  = @(thing) (~isempty(thing) & ~sum(isnan(thing)));
    p       = sys.p;
    q       = sys.q;
    % convert x+/-e into precision for p(x) (precision = 1/variance)
    w       = @(x,e) abs( p(1 + e./x) ).^(-2);

    nv      = size(sys.K,2);
    yobs    = nan(nD,nv);
    wobs    = nan(nD,nv);

    if (~isfield(obs,'tp'))
        if opt.printmes ~= 0
            error('Need to provide temperature and pressure measurement.')
        end
    end
    
    nTP = length(obs(1).tp); 
    for i = 1:nD % loop over all the stations
        % make sure all the required fields in the obs struct exist
        if (~isfield(obs(i), 'sal'))
            if opt.printmes ~= 0
                error('Need to provide salinity measurement.');
            end
        else
            yobs(i,sys.isal) = obs(i).sal;
        end
        if (~isfield(obs(i), 'esal'))
            obs(i).esal         = 0.002; % std = 0.002 PSU
            wobs(i,sys.isal)    = (obs(i).esal)^(-2);
            if opt.printmes ~= 0
                fprintf('Warning: Assuming salinity uncertainty is 0.002 PSU \n');
            end
        else
            wobs(i,sys.isal)    = (obs(i).esal)^(-2); % std e -> w
        end
        if (~isfield(obs(i),'TC')) || (~isgood(obs(i).TC))
            obs(i).TC        = nan; %[]
            yobs(i,sys.ipTC) = nan;
        else
            yobs(i,sys.ipTC) = p((obs(i).TC)*1e-6); % convt to mol/kg
        end
        if (~isfield(obs(i),'eTC')) || (~isgood(obs(i).eTC))
            obs(i).eTC       = nan;
            wobs(i,sys.ipTC) = nan;
        else
            wobs(i,sys.ipTC) = w(obs(i).TC,obs(i).eTC); % std e -> w
        end
        if(~isfield(obs(i),'TA'))  || (~isgood(obs(i).TA))
            obs(i).TA        = nan; %[]
            yobs(i,sys.ipTA) = nan;
        else
            yobs(i,sys.ipTA) = p((obs(i).TA)*1e-6); % convt to mol/kg
        end
        if (~isfield(obs(i),'eTA'))  || (~isgood(obs(i).eTA))
            obs(i).eTA       = nan;
            wobs(i,sys.ipTA) = nan;
        else
            wobs(i,sys.ipTA) = w(obs(i).TA,obs(i).eTA); % std e -> w
        end
        % calculate totals that are a function of salinity
        [pT,~,~,epT]  =  calc_pTOT(opt,obs(i).sal);
        pTB  = pT(1); epTB  = epT(1);
        pTS  = pT(2); epTS  = epT(2);
        pTF  = pT(3); epTF  = epT(3);
        pTCa = pT(4); epTCa = epT(4); % (see Ref's within calc_pTOT)
        % total borate
        if (~isfield(obs(i),'TB')) || (~isgood(obs(i).TB))
            obs(i).TB = nan;
            if opt.turnoff.TB == 1 % do not input TB
                yobs(i,sys.ipTB) = nan; 
            else % do input TB
                yobs(i,sys.ipTB) = pTB;
            end
        else
            if ((obs(i).TB) == 0)
                obs(i).TB   = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTB) = p(obs(i).TB*1e-6); % convt µmol/kg to mol/kg
        end
        if (~isfield(obs(i), 'eTB'))  || (~isgood(obs(i).eTB))
            obs(i).eTB       = nan;
            if opt.turnoff.TB == 1 % do not input eTB
                wobs(i,sys.ipTB) = nan;
            else % do input eTB
                wobs(i,sys.ipTB) = (epTB)^(-2); % convert to precision
            end
        else
            wobs(i,sys.ipTB) = w(obs(i).TB,obs(i).eTB); % mol/kg
        end
        
        % total sulfate
        if (~isfield(obs(i), 'TS'))  || (~isgood(obs(i).TS))
            obs(i).TS        = nan;
            yobs(i,sys.ipTS) = pTS;
        else
            if ((obs(i).TS) == 0)
                obs(i).TS = 1e-3; % µmol/kg reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTS)     = p(obs(i).TS*1e-6); % mol/kg
        end
        if (~isfield(obs(i), 'eTS'))  || (~isgood(obs(i).eTS))
            obs(i).eTS       = nan ;
            wobs(i,sys.ipTS) = (epTS)^(-2); % convert to precision
        else
            TS = (0.14/96.062)*obs(i).sal/1.80655;
            wobs(i,sys.ipTS)    = w(TS,obs(i).eTS);
            % wobs(i,sys.ipTS)    = w(obs(i).TS,obs(i).eTS);
        end
        % total fluoride
        if (~isfield(obs(i), 'TF'))  || (~isgood(obs(i).TF))
            obs(i).TF        = nan; 
            yobs(i,sys.ipTF) = pTF;
        else
            if ((obs(i).TF) == 0)
                obs(i).TF       = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTF)    = p(obs(i).TF*1e-6); % convt µmol/kg to mol/kg
        end
        if (~isfield(obs(i), 'eTF'))  || (~isgood(obs(i).eTF))
            obs(i).eTF       = nan;
            wobs(i,sys.ipTF) = (epTF)^(-2);
        else
            TF = (6.7e-5/18.998)*obs(i).sal/1.80655;
            wobs(i,sys.ipTF)    = w(TF,obs(i).eTF);
            % wobs(i,sys.ipTF)    = w(obs(i).TF,obs(i).eTF);
        end        
        
        % total phosphate
        if (~isfield(obs(i), 'TP'))  || (~isgood(obs(i).TP))
            obs(i).TP           = nan;
            yobs(i,sys.ipTP)    = p(1e-9); % resest minimum (mol/kg)
        else
            if ((obs(i).TP) == 0) % zero po4 is very unlikely and breaks the code
                obs(i).TP       = 1e-3; % umol/kg-SW, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTP)     = p(obs(i).TP*1e-6); % convt µmol/kg to mol/kg
        end
        if (~isfield(obs(i), 'eTP'))  || (~isgood(obs(i).eTP))
            wobs(i,sys.ipTP)    = w(1e-3,1e-3); % mol/kg
            obs(i).eTP          = nan;
        else
            if ((obs(i).eTP) == 0)
                obs(i).eTP      = 1e-3; % umol/kg, reset minimum if zero
            end
            wobs(i,sys.ipTP)    = w(obs(i).TP,obs(i).eTP); % mol/kg
        end
        
        % total silicate
        if (~isfield(obs(i), 'TSi'))  || (~isgood(obs(i).TSi))
            obs(i).TSi          = nan;
            yobs(i,sys.ipTSi)   = p(1e-9); % reset minimum (mol/kg)
        else
            if ((obs(i).TSi) == 0) % zero silicate very unlikely and breaks code
                obs(i).TSi      = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTSi)   = p(obs(i).TSi*1e-6);
        end
        if (~isfield(obs(i), 'eTSi'))  || (~isgood(obs(i).eTSi))
            wobs(i,sys.ipTSi)   = w(1e-3,1e-3); % mol/kg
            obs(i).eTSi         = nan;
            if (isgood(obs(i).TSi)) && opt.printmes ~= 0
                fprintf('Warning, no obs.eTSi input with obs.eTSi. Assuming 1 nanomolar.\n' )
            end
        else
            if ((obs(i).eTSi) == 0)
                obs(i).eTSi     = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            wobs(i,sys.ipTSi)   = w(obs(i).TSi,obs(i).eTSi); %  mol/kg
        end
        % total amonia
        if (~isfield(obs(i), 'TNH4'))  || (~isgood(obs(i).TNH4))
            obs(i).TNH4         = nan;
            yobs(i,sys.ipTNH4)  = p(1e-9); % reset minimum (mol/kg)
        else
            if ((obs(i).TNH4) == 0)
                obs(i).TNH4     = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTNH4)  = p(obs(i).TNH4*1e-6);
        end
        if (~isfield(obs(i), 'eTNH4'))  || (~isgood(obs(i).eTNH4))
            eTNH4               = 5e-4; % µmol/kg
            wobs(i,sys.ipTNH4)  = w(1e-3,eTNH4); % mol/kg
            obs(i).eTNH4        = nan;
            if (isgood(obs(i).TNH4)) && opt.printmes ~= 0
                fprintf('Warning, no obs.eTNH4 input with obs.eTNH4. Assuming 5e-4 umol/kg.\n' )
            end
        else
            wobs(i,sys.ipTNH4)  = w(obs(i).TNH4,obs(i).eTNH4);
        end
        % total sulfide
        if (~isfield(obs(i), 'TH2S'))  || (~isgood(obs(i).TH2S))
            obs(i).TH2S         = nan; 
            yobs(i,sys.ipTH2S)  = p(1e-9); % reset minimum (mol/kg)
        else
            if ((obs(i).TH2S) == 0)
                obs(i).TH2S     = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTH2S)  = p(obs(i).TH2S*1e-6);
        end
        if (~isfield(obs(i), 'eTH2S'))  || (~isgood(obs(i).eTH2S))
            eTH2S               = 5e-4; % µmol/kg
            wobs(i,sys.ipTH2S)  = w(1e-3,eTH2S); % mol/kg
            obs(i).eTH2S        = nan;
            if (isgood(obs(i).TH2S)) && opt.printmes ~= 0
                fprintf(' Warning, no obs.eTH2S input with obs.eTNH4. Assuming 5e-4 umol/kg.\n' )
            end
        else
            wobs(i,sys.ipTH2S)  = w(obs(i).TH2S,obs(i).eTH2S);
        end
            
        % total calcium
        if (~isfield(obs(i), 'TCa'))  || (~isgood(obs(i).TCa))
            obs(i).TCa        = nan;
            yobs(i,sys.ipTCa) = pTCa;
        else
            if ((obs(i).TCa) == 0)
                obs(i).TCa      = 1e-9; % mol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTCa)   = p(obs(i).TCa); % assume user input of mol/kg
        end
        if (~isfield(obs(i), 'eTCa'))  || (~isgood(obs(i).eTCa))
            obs(i).eTCa       = nan;
            wobs(i,sys.ipTCa) = (epTCa)^(-2);
            if (isgood(obs(i).TCa)) && opt.printmes ~= 0
                fprintf(' Warning, no obs.eTCa input with obs.eTCa. Assuming 6e-5 mol/kg.\n' )
            end
        else
            wobs(i,sys.ipTCa)   = w(obs(i).TCa,obs(i).eTCa);
        end
        
        if opt.pKalpha == 1
            % Organic Alkalinity Alpha
            if (~isfield(obs(i), 'TAlpha')) || (~isgood(obs(i).TAlpha))
                obs(i).TAlpha           = nan; 
                yobs(i,sys.ipTAlpha)    = p(5e-6); % guess 5 umol/kg
            else
                if ((obs(i).TAlpha) == 0)
                    obs(i).TAlpha       = 1e-9;
                end
                yobs(i,sys.ipTAlpha)    = p(obs(i).TAlpha*1e-6);
            end
            if (~isfield(obs(i), 'eTAlpha')) || (~isgood(obs(i).eTAlpha))
                obs(i).eTAlpha          = nan; 
                wobs(i,sys.ipTAlpha)    = w(5,2); % 5 ± 2 umol/kg
            else
                wobs(i,sys.ipTAlpha)    = w(obs(i).TAlpha,obs(i).eTAlpha);
            end
        end

        if opt.pKbeta == 1
            % Organic Alkalinity Beta
            if (~isfield(obs(i),'TBeta')) || (~isgood(obs(i).TBeta))
                obs(i).TBeta        = nan;
                yobs(i,sys.ipTBeta) = p(5e-6);
            else
                if ((obs(i).TBeta) == 0)
                    obs(i).TBeta    = 1e-9;
                end
                yobs(i,sys.ipTBeta) = p(obs(i).TBeta*1e-6);
            end
            if (~isfield(obs(i),'eTBeta')) || (~isgood(obs(i).eTBeta))
                obs(i).eTBeta       = nan;
                wobs(i,sys.ipTBeta) = w(5,2);
            else
                wobs(i,sys.ipTBeta) = w(obs(i).TBeta,obs(i).eTBeta);
            end
        end

        for j = 1:nTP % loop over (T,P) systems
            yobs(i,sys.tp(j).iT)   = obs(i).tp(j).T;
            yobs(i,sys.tp(j).iP)   = obs(i).tp(j).P;
            wobs(i,sys.tp(j).iT)   = (obs(i).tp(j).eT)^(-2);
            wobs(i,sys.tp(j).iP)   = (obs(i).tp(j).eP)^(-2);
            [pK,~,epK] = calc_pK(opt,obs(i).tp(j).T,obs(i).sal,obs(i).tp(j).P); % T, S, P

            pK0   = pK(1);      pK1  = pK(2);     pK2   = pK(3);  
            pKb   = pK(4);      pKw  = pK(5);     pKs   = pK(6);  
            pKf   = pK(7);      pKp1 = pK(8);     pKp2  = pK(9);  
            pKp3  = pK(10);     pKsi = pK(11);    pKnh4 = pK(12);
            pKh2s = pK(13);     pp2f = pK(14);    pKar  = pK(15); 
            pKca  = pK(16);     pfH  = pK(17);
            
            epK0   = epK(1);    epK1  = epK(2);   epK2   = epK(3);  
            epKb   = epK(4);    epKw  = epK(5);   epKs   = epK(6);  
            epKf   = epK(7);    epKp1 = epK(8);   epKp2  = epK(9);  
            epKp3  = epK(10);   epKsi = epK(11);  epKnh4 = epK(12);
            epKh2s = epK(13);   epp2f = epK(14);  epKar  = epK(15); 
            epKca  = epK(16);   epfH  = epK(17);
            
            % add "observations" for the equilibrium constants
            % and transfer from obs struct to yobs and wobs
            
            % co2 solubility and fugacity
            if (~isfield(obs(i).tp(j),'pK0')) || (~isgood(obs(i).tp(j).pK0))
                obs(i).tp(j).pK0        = nan;
                yobs(i,sys.tp(j).ipK0)  = pK0;
            else
                yobs(i,sys.tp(j).ipK0)  = obs(i).tp(j).pK0;
            end
            if (~isfield(obs(i).tp(j),'epK0')) || (~isgood(obs(i).tp(j).epK0))
                obs(i).tp(j).epK0       = nan;
                wobs(i,sys.tp(j).ipK0)  = (epK0)^(-2);
            else
                wobs(i,sys.tp(j).ipK0)  = (obs(i).tp(j).pK0)^(-2);
            end
            if (~isfield(obs(i).tp(j),'co2st')) || (~isgood(obs(i).tp(j).co2st))
                obs(i).tp(j).co2st          = nan;
                yobs(i,sys.tp(j).ipco2st)   = nan;
            else
                yobs(i,sys.tp(j).ipco2st)   = p(obs(i).tp(j).co2st*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'eco2st')) || (~isgood(obs(i).tp(j).eco2st))
                obs(i).tp(j).eco2st         = nan;
                wobs(i,sys.tp(j).ipco2st)   = nan;
            else
                wobs(i,sys.tp(j).ipco2st) = w(obs(i).tp(j).co2st,obs(i).tp(j).eco2st);
            end
            if (~isfield(obs(i).tp(j),'fco2')) || (~isgood(obs(i).tp(j).fco2))
                obs(i).tp(j).fco2           = nan;
                yobs(i,sys.tp(j).ipfco2)    = nan;
            else
                yobs(i,sys.tp(j).ipfco2)    = p(obs(i).tp(j).fco2*1e-6); % convt µatm to atm
            end
            if (~isfield(obs(i).tp(j),'efco2')) || (~isgood(obs(i).tp(j).efco2))
                obs(i).tp(j).efco2          = nan;
                wobs(i,sys.tp(j).ipfco2)    = nan;
            else
                wobs(i,sys.tp(j).ifco2) = w(obs(i).tp(j).fco2, obs(i).tp(j).efco2);
            end
            if (~isfield(obs(i).tp(j),'pp2f')) || (~isgood(obs(i).tp(j).pp2f))
                obs(i).tp(j).pp2f           = nan;
                yobs(i,sys.tp(j).ipp2f)     = pp2f;
            else
                yobs(i,sys.tp(j).ipp2f)     = obs(i).tp(j).pp2f;
            end
            if (~isfield(obs(i).tp(j),'epp2f')) || (~isgood(obs(i).tp(j).epp2f))
                obs(i).tp(j).epp2f         = nan;
                wobs(i,sys.tp(j).ipp2f)    = (epp2f)^(-2);
            else
                wobs(i,sys.tp(j).ipp2f) = (obs(i).tp(j).epp2f)^(-2);
            end
            if (~isfield(obs(i).tp(j),'pco2')) || (~isgood(obs(i).tp(j).pco2))
                obs(i).tp(j).pco2          = nan;
                yobs(i,sys.tp(j).ippco2)   = nan;
            else
                yobs(i,sys.tp(j).ippco2)   = p(obs(i).tp(j).pco2*1e-6); % convt µatm to atm
            end
            if (~isfield(obs(i).tp(j),'epco2')) || (~isgood(obs(i).tp(j).epco2))
                obs(i).tp(j).epco2         = nan;
                wobs(i,sys.tp(j).ippco2)   = nan;
            else
                wobs(i,sys.tp(j).ippco2) = w(obs(i).tp(j).pco2,obs(i).tp(j).epco2);
            end

            % carbonate system
            if (~isfield(obs(i).tp(j),'pK1')) || (~isgood(obs(i).tp(j).pK1))
                obs(i).tp(j).pK1       = nan;
                if opt.turnoff.pK1 == 1 % do not input pK1
                    if nTP > 1
                        error(['opt.turnoff.pK1 only works with one tp, not 2+. ' ...
                                'Must input at tp(1) only.']);
                    else
                        yobs(i,sys.tp(j).ipK1) = nan;
                    end
                else % do input pK1
                    yobs(i,sys.tp(j).ipK1) = pK1;
                end
            else
                yobs(i,sys.tp(j).ipK1) = obs(i).tp(j).pK1;
            end
            if (~isfield(obs(i).tp(j),'epK1')) || (~isgood(obs(i).tp(j).epK1))
                obs(i).tp(j).epK1      = nan;
                if opt.turnoff.pK1 == 1
                    wobs(i,sys.tp(j).ipK1) = nan;
                else
                    wobs(i,sys.tp(j).ipK1) = (epK1)^(-2);
                end
            else
                wobs(i,sys.tp(j).ipK1) = (obs(i).tp(j).epK1)^(-2);
            end
            if (~isfield(obs(i).tp(j),'pK2')) || (~isgood(obs(i).tp(j).pK2))
                obs(i).tp(j).pK2       = nan;
                if opt.turnoff.pK2 == 1
                    if nTP > 1
                        error(['opt.turnoff.pK2 only works with one tp, not 2+. ' ...
                                'Must input at tp(1) only.']);
                    else
                        yobs(i,sys.tp(j).ipK2) = nan;
                    end
                else
                    yobs(i,sys.tp(j).ipK2) = pK2;
                end
            else
                yobs(i,sys.tp(j).ipK2) = obs(i).tp(j).pK2;
            end
            if (~isfield(obs(i).tp(j),'epK2')) || (~isgood(obs(i).tp(j).epK2))
                obs(i).tp(j).epK2      = nan;
                if opt.turnoff.pK2 == 1
                    wobs(i,sys.tp(j).ipK2) = nan;
                else
                    wobs(i,sys.tp(j).ipK2) = (epK2)^(-2);
                end
            else
                wobs(i,sys.tp(j).ipK2) = (obs(i).tp(j).epK2)^(-2);
            end
            if (~isfield(obs(i).tp(j),'hco3')) || (~isgood(obs(i).tp(j).hco3))
                obs(i).tp(j).hco3          = nan;
                yobs(i,sys.tp(j).iphco3)   = nan;   
            else
                yobs(i,sys.tp(j).iphco3) = p(obs(i).tp(j).hco3*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'ehco3')) || (~isgood(obs(i).tp(j).ehco3))
                obs(i).tp(j).ehco3         = nan;
                wobs(i,sys.tp(j).iphco3)   = nan;
            else
                wobs(i,sys.tp(j).iphco3) = w(obs(i).tp(j).hco3,obs(i).tp(j).ehco3);
            end
            if (~isfield(obs(i).tp(j),'co3')) || (~isgood(obs(i).tp(j).co3))
                obs(i).tp(j).co3           = nan;
                yobs(i,sys.tp(j).ipco3)    = nan;
            else
                yobs(i,sys.tp(j).ipco3) = p(obs(i).tp(j).co3*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'eco3')) || (~isgood(obs(i).tp(j).eco3))
                obs(i).tp(j).eco3          = nan;
                wobs(i,sys.tp(j).ipco3)    = nan;
            else
                wobs(i,sys.tp(j).ipco3) = w(obs(i).tp(j).co3,obs(i).tp(j).eco3);
            end
            if (~isfield(obs(i).tp(j),'ph')) || (~isgood(obs(i).tp(j).ph))
                obs(i).tp(j).ph        = nan;
                yobs(i,sys.tp(j).iph)  = nan;
            else
                yobs(i,sys.tp(j).iph)  = obs(i).tp(j).ph;
            end
            if (~isfield(obs(i).tp(j),'eph')) || (~isgood(obs(i).tp(j).eph))
                obs(i).tp(j).eph       = nan;
                wobs(i,sys.tp(j).iph)  = nan;
            else
                wobs(i,sys.tp(j).iph)  = (obs(i).tp(j).eph).^(-2);
            end
            if (~isfield(obs(i).tp(j),'ph_tot')) || (~isgood(obs(i).tp(j).ph_tot))
                obs(i).tp(j).ph_tot       = nan;
                yobs(i,sys.tp(j).iph_tot) = nan;
            else
                yobs(i,sys.tp(j).iph_tot) = obs(i).tp(j).ph_tot;
            end
            if (~isfield(obs(i).tp(j),'eph_tot')) || (~isgood(obs(i).tp(j).eph_tot))
                obs(i).tp(j).eph_tot      = nan;
                wobs(i,sys.tp(j).iph_tot) = nan;
            else
                wobs(i,sys.tp(j).iph_tot) = obs(i).tp(j).eph_tot.^(-2);
            end
            if (~isfield(obs(i).tp(j),'ph_free')) || (~isgood(obs(i).tp(j).ph_free))
                obs(i).tp(j).ph_free       = nan;
                yobs(i,sys.tp(j).iph_free) = nan;
            else
                yobs(i,sys.tp(j).iph_free) = obs(i).tp(j).ph_free ;
            end
            if (~isfield(obs(i).tp(j),'eph_free')) || (~isgood(obs(i).tp(j).eph_free))
                obs(i).tp(j).eph_free      = nan;
                wobs(i,sys.tp(j).iph_free) = nan;
            else
                wobs(i,sys.tp(j).iph_free) = obs(i).tp(j).eph_free.^(-2);
            end
            if (~isfield(obs(i).tp(j),'ph_sws')) || (~isgood(obs(i).tp(j).ph_sws))
                obs(i).tp(j).ph_sws       = nan;
                yobs(i,sys.tp(j).iph_sws) = nan;
            else
                yobs(i,sys.tp(j).iph_sws) = obs(i).tp(j).ph_sws;
            end
            if (~isfield(obs(i).tp(j),'eph_sws')) || (~isgood(obs(i).tp(j).eph_sws))
                obs(i).tp(j).eph_sws      = nan;
                wobs(i,sys.tp(j).iph_sws) = nan;
            else
                wobs(i,sys.tp(j).iph_sws) = obs(i).tp(j).eph_sws.^(-2);
            end
            if (~isfield(obs(i).tp(j),'ph_nbs')) || (~isgood(obs(i).tp(j).ph_nbs))
                obs(i).tp(j).ph_nbs       = nan;
                yobs(i,sys.tp(j).iph_nbs) = nan;
            else
                yobs(i,sys.tp(j).iph_nbs) = obs(i).tp(j).ph_nbs;
            end
            if (~isfield(obs(i).tp(j),'eph_nbs')) || (~isgood(obs(i).tp(j).eph_nbs))
                obs(i).tp(j).eph_nbs      = nan;
                wobs(i,sys.tp(j).iph_nbs) = nan;
            else
                wobs(i,sys.tp(j).iph_nbs) = obs(i).tp(j).eph_nbs.^(-2);
            end
            if (~isfield(obs(i).tp(j),'pfH')) || (~isgood(obs(i).tp(j).pfH))
                obs(i).tp(j).pfH           = nan;
                yobs(i,sys.tp(j).ipfH)     = pfH;
            else
                yobs(i,sys.tp(j).ipfH)     = obs(i).tp(j).pfH;
            end
            if (~isfield(obs(i).tp(j),'epfH')) || (~isgood(obs(i).tp(j).epfH))
                obs(i).tp(j).epfH          = nan;
                wobs(i,sys.tp(j).ipfH)     = (epfH).^(-2);
            else
                wobs(i,sys.tp(j).ipfH)     = (obs(i).tp(j).epfH).^(-2) ;
            end
            
            % water dissociation
            if (~isfield(obs(i).tp(j),'pKw')) || (~isgood(obs(i).tp(j).pKw))
                obs(i).tp(j).pKw       = nan;
                yobs(i,sys.tp(j).ipKw) = pKw;
            else
                yobs(i,sys.tp(j).ipKw) = obs(i).tp(j).pKw;
            end
            if (~isfield(obs(i).tp(j),'epKw')) || (~isgood(obs(i).tp(j).epKw))
                obs(i).tp(j).epKw      = nan;
                wobs(i,sys.tp(j).ipKw) = (epKw).^(-2);
            else
                wobs(i,sys.tp(j).ipKw) = (obs(i).tp(j).epKw).^(-2);
            end
            if (~isfield(obs(i).tp(j),'oh')) || (~isgood(obs(i).tp(j).oh))
                obs(i).tp(j).oh        = nan;
                yobs(i,sys.tp(j).ipoh) = nan;
            else
                yobs(i,sys.tp(j).ipoh) = p(obs(i).tp(j).oh*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'eoh')) || (~isgood(obs(i).tp(j).eoh))
                obs(i).tp(j).eoh       = nan;
                wobs(i,sys.tp(j).ipoh) = nan;
            else
                wobs(i,sys.tp(j).ipoh) = w(obs(i).tp(j).oh,obs(i).tp(j).eoh);
            end

            % borate system 
            if (~isfield(obs(i).tp(j),'pKb')) || (~isgood(obs(i).tp(j).pKb))
                obs(i).tp(j).pKb       = nan;
                yobs(i,sys.tp(j).ipKb) = pKb; 
            else
                yobs(i,sys.tp(j).ipKb) = obs(i).tp(j).pKb;
            end
            if (~isfield(obs(i).tp(j),'epKb')) || (~isgood(obs(i).tp(j).epKb))
                obs(i).tp(j).epKb      = nan;
                wobs(i,sys.tp(j).ipKb) = (epKb).^(-2);
            else
                wobs(i,sys.tp(j).ipKb) = (obs(i).tp(j).epKb).^(-2);
            end
            if (~isfield(obs(i).tp(j),'boh3')) || (~isgood(obs(i).tp(j).boh3))
                obs(i).tp(j).boh3          = nan;
                yobs(i,sys.tp(j).ipboh3)   = nan;
            else
                yobs(i,sys.tp(j).ipboh3) = p(obs(i).tp(j).boh3*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'eboh3')) || (~isgood(obs(i).tp(j).eboh3))
                obs(i).tp(j).eboh3         = nan;
                wobs(i,sys.tp(j).ipboh3)   = nan;
            else
                wobs(i,sys.tp(j).ipboh3) = w(obs(i).tp(j).boh3,obs(i).tp(j).eboh3);
            end
            if (~isfield(obs(i).tp(j),'boh4')) || (~isgood(obs(i).tp(j).boh4))
                obs(i).tp(j).boh4          = nan;
                yobs(i,sys.tp(j).ipboh4)   = nan;
            else
                yobs(i,sys.tp(j).ipboh4) = p(obs(i).tp(j).boh4*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'eboh4')) || (~isgood(obs(i).tp(j).eboh4))
                obs(i).tp(j).eboh4         = nan;
                wobs(i,sys.tp(j).ipboh4)   = nan;
            else
                wobs(i,sys.tp(j).ipboh4) = w(obs(i).tp(j).boh4,obs(i).tp(j).eboh4);
            end
            
            % sulfate system
            if (~isfield(obs(i).tp(j),'pKs')) || (~isgood(obs(i).tp(j).pKs))
                obs(i).tp(j).pKs       = nan;
                yobs(i,sys.tp(j).ipKs) = pKs;
            else
                yobs(i,sys.tp(j).ipKs) = obs(i).tp(j).pKs;
            end
            if (~isfield(obs(i).tp(j),'epKs')) || (~isgood(obs(i).tp(j).epKs))
                obs(i).tp(j).epKs      = nan;
                wobs(i,sys.tp(j).ipKs) = (epKs).^(-2);
            else
                wobs(i,sys.tp(j).ipKs) = (obs(i).tp(j).epKs).^(-2);
            end
            if (~isfield(obs(i).tp(j),'hso4')) || (~isgood(obs(i).tp(j).hso4))
                obs(i).tp(j).hso4          = nan;
                yobs(i,sys.tp(j).iphso4)   = nan;
            else
                yobs(i,sys.tp(j).iphso4) = p(obs(i).tp(j).hso4*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'ehso4')) || (~isgood(obs(i).tp(j).ehso4))
                obs(i).tp(j).ehso4         = nan;
                wobs(i,sys.tp(j).iphso4)   = nan;
            else
                wobs(i,sys.tp(j).iphso4) = w(obs(i).tp(j).hso4,obs(i).tp(j).ehso4);
            end
            if (~isfield(obs(i).tp(j),'so4')) || (~isgood(obs(i).tp(j).so4))
                obs(i).tp(j).so4           = nan;
                yobs(i,sys.tp(j).ipso4)    = nan;
            else
                yobs(i,sys.tp(j).ipso4) = p(obs(i).tp(j).so4*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'eso4')) || (~isgood(obs(i).tp(j).eso4))
                obs(i).tp(j).eso4          = nan;
                wobs(i,sys.tp(j).ipso4)    = nan;
            else
                wobs(i,sys.tp(j).ipso4) = w(obs(i).tp(j).so4,obs(i).tp(j).eso4);
            end

            % fluoride system
            if (~isfield(obs(i).tp(j),'pKf')) || (~isgood(obs(i).tp(j).pKf))
                obs(i).tp(j).pKf       = nan;
                yobs(i,sys.tp(j).ipKf) = pKf;
            else
                yobs(i,sys.tp(j).ipKf) = obs(i).tp(j).pKf;
            end
            if (~isfield(obs(i).tp(j),'epKf')) || (~isgood(obs(i).tp(j).epKf))
                obs(i).tp(j).epKf      = nan;
                wobs(i,sys.tp(j).ipKf) = (epKf).^(-2);
            else
                wobs(i,sys.tp(j).ipKf) = (obs(i).tp(j).epKf).^(-2);
            end
            if (~isfield(obs(i).tp(j),'HF')) || (~isgood(obs(i).tp(j).HF))
                obs(i).tp(j).HF        = nan;
                yobs(i,sys.tp(j).ipHF) = nan;
            else
                yobs(i,sys.tp(j).ipHF) = p(obs(i).tp(j).HF*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'eHF')) || (~isgood(obs(i).tp(j).eHF))
                obs(i).tp(j).eHF       = nan;
                wobs(i,sys.tp(j).ipHF) = nan;
            else
                wobs(i,sys.tp(j).ipHF) = w(obs(i).tp(j).HF,obs(i).tp(j).eHF);
            end
            if (~isfield(obs(i).tp(j),'F')) || (~isgood(obs(i).tp(j).F))
                obs(i).tp(j).F         = nan;
                yobs(i,sys.tp(j).ipF)  = nan;
            else
                yobs(i,sys.tp(j).ipF) = p(obs(i).tp(j).F*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'eF')) || (~isgood(obs(i).tp(j).eF))
                obs(i).tp(j).eF        = nan;
                wobs(i,sys.tp(j).ipF)  = nan;
            else
                wobs(i,sys.tp(j).ipF) = w(obs(i).tp(j).F,obs(i).tp(j).eF);
            end
            
            % phosphate system
            if (~isfield(obs(i).tp(j),'pKp1')) || (~isgood(obs(i).tp(j).pKp1))
                obs(i).tp(j).pKp1          = nan;
                yobs(i,sys.tp(j).ipKp1)    = pKp1;
            else
               yobs(i,sys.tp(j).ipKp1)    = obs(i).tp(j).pKp1; 
            end
            if (~isfield(obs(i).tp(j),'epKp1')) || (~isgood(obs(i).tp(j).epKp1))
                obs(i).tp(j).epKp1         = nan;
                wobs(i,sys.tp(j).ipKp1)    = (epKp1).^(-2);
            else
                wobs(i,sys.tp(j).ipKp1) = (obs(i).tp(j).epKp1).^(-2);
            end
            if (~isfield(obs(i).tp(j),'pKp2')) || (~isgood(obs(i).tp(j).pKp2))
               obs(i).tp(j).pKp2          = nan;
                yobs(i,sys.tp(j).ipKp2)    = pKp2; 
            else
                yobs(i,sys.tp(j).ipKp2)    = obs(i).tp(j).pKp2;
            end
            if (~isfield(obs(i).tp(j),'epKp2')) || (~isgood(obs(i).tp(j).epKp2))
                obs(i).tp(j).epKp2         = nan;
                wobs(i,sys.tp(j).ipKp2)    = (epKp2).^(-2);
            else
                wobs(i,sys.tp(j).ipKp2) = (obs(i).tp(j).epKp2).^(-2);
            end
            if (~isfield(obs(i).tp(j),'pKp3')) || (~isgood(obs(i).tp(j).pKp3))
                obs(i).tp(j).pKp3          = nan;
                yobs(i,sys.tp(j).ipKp3)    = pKp3;
            else
                yobs(i,sys.tp(j).ipKp3)    = obs(i).tp(j).pKp3;
            end
            if (~isfield(obs(i).tp(j),'epKp3')) || (~isgood(obs(i).tp(j).epKp3))
                obs(i).tp(j).epKp3         = nan;
                wobs(i,sys.tp(j).ipKp3)    = (epKp3).^(-2);
            else
                wobs(i,sys.tp(j).ipKp3) = (obs(i).tp(j).epKp3).^(-2);
            end        
            if (~isfield(obs(i).tp(j),'h3po4')) || (~isgood(obs(i).tp(j).h3po4))
                obs(i).tp(j).h3po4         = nan;
                yobs(i,sys.tp(j).iph3po4)  = nan;
            else
                yobs(i,sys.tp(j).iph3po4) = p(obs(i).tp(j).h3po4*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'eh3po4')) || (~isgood(obs(i).tp(j).eh3po4))
                obs(i).tp(j).eh3po4        = nan;
                wobs(i,sys.tp(j).iph3po4)  = nan;
            else
                wobs(i,sys.tp(j).iph3po4) = w(obs(i).tp(j).h3po4,obs(i).tp(j).eh3po4);
            end
            if (~isfield(obs(i).tp(j),'h2po4')) || (~isgood(obs(i).tp(j).h2po4))
                obs(i).tp(j).h2po4         = nan;
                yobs(i,sys.tp(j).iph2po4)  = nan;
            else
                yobs(i,sys.tp(j).iph2po4) = p(obs(i).tp(j).h2po4*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'eh2po4')) || (~isgood(obs(i).tp(j).eh2po4))
                obs(i).tp(j).eh2po4        = nan;
                wobs(i,sys.tp(j).iph2po4)  = nan;
            else
                wobs(i,sys.tp(j).iph2po4) = w(obs(i).tp(j).h2po4,obs(i).tp(j).eh2po4);
            end
            if (~isfield(obs(i).tp(j),'hpo4')) || (~isgood(obs(i).tp(j).hpo4))
                obs(i).tp(j).hpo4          = nan;
                yobs(i,sys.tp(j).iphpo4)   = nan;
            else
               yobs(i,sys.tp(j).iphpo4) = p(obs(i).tp(j).hpo4*1e-6); % convt µmol/kg to mol/kg 
            end
            if (~isfield(obs(i).tp(j),'ehpo4')) || (~isgood(obs(i).tp(j).ehpo4))
                obs(i).tp(j).ehpo4         = nan;
                wobs(i,sys.tp(j).iphpo4)   = nan;
            else
                wobs(i,sys.tp(j).iphpo4) = w(obs(i).tp(j).hpo4,obs(i).tp(j).ehpo4);
            end
            if (~isfield(obs(i).tp(j),'po4')) || (~isgood(obs(i).tp(j).po4))
                obs(i).tp(j).po4           = nan;
                yobs(i,sys.tp(j).ippo4)    = nan;
            else
                yobs(i,sys.tp(j).ippo4) = p(obs(i).tp(j).po4*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'epo4')) || (~isgood(obs(i).tp(j).epo4))
                obs(i).tp(j).epo4          = nan;
                wobs(i,sys.tp(j).ippo4)    = nan;
            else
                wobs(i,sys.tp(j).ippo4) = w(obs(i).tp(j).po4,obs(i).tp(j).epo4);
            end

            % silicate system
            if (~isfield(obs(i).tp(j),'pKsi')) || (~isgood(obs(i).tp(j).pKsi))
                obs(i).tp(j).pKsi          = nan;
                yobs(i,sys.tp(j).ipKsi)    = pKsi;
            else
                yobs(i,sys.tp(j).ipKsi)    = obs(i).tp(j).pKsi;
            end
            if (~isfield(obs(i).tp(j),'epKsi')) || (~isgood(obs(i).tp(j).epKsi))
                obs(i).tp(j).epKsi         = nan;
                wobs(i,sys.tp(j).ipKsi)    = (epKsi).^(-2); 
            else
                wobs(i,sys.tp(j).ipKsi) = (obs(i).tp(j).epKsi).^(-2);
            end
            if (~isfield(obs(i).tp(j),'siooh3')) || (~isgood(obs(i).tp(j).siooh3))
                obs(i).tp(j).siooh3        = nan;
                yobs(i,sys.tp(j).ipsiooh3) = nan;
            else
                yobs(i,sys.tp(j).ipsiooh3) = p(obs(i).tp(j).siooh3*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'esiooh3')) || (~isgood(obs(i).tp(j).esiooh3))
                obs(i).tp(j).esiooh3       = nan;
                wobs(i,sys.tp(j).ipsiooh3) = nan;
            else
                wobs(i,sys.tp(j).ipsiooh3) = w(obs(i).tp(j).siooh3,obs(i).tp(j).esiooh3);
            end
            if (~isfield(obs(i).tp(j),'sioh4')) || (~isgood(obs(i).tp(j).sioh4))
                obs(i).tp(j).sioh4         = nan;
                yobs(i,sys.tp(j).ipsioh4)  = nan;
            else
                yobs(i,sys.tp(j).ipsioh4) = p(obs(i).tp(j).sioh4*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'esioh4')) || (~isgood(obs(i).tp(j).esioh4))
                obs(i).tp(j).esioh4        = nan;
                wobs(i,sys.tp(j).ipsioh4)  = nan;
            else
                wobs(i,sys.tp(j).ipsioh4) = w(obs(i).tp(j).sioh4,obs(i).tp(j).esioh4);
            end

            % ammonia system
            if (~isfield(obs(i).tp(j),'pKnh4')) || (~isgood(obs(i).tp(j).pKnh4))
                obs(i).tp(j).pKnh4         = nan;
                yobs(i,sys.tp(j).ipKnh4)   = pKnh4;
            else
                yobs(i,sys.tp(j).ipKnh4)   = obs(i).tp(j).pKnh4;
            end
            if (~isfield(obs(i).tp(j),'epKnh4')) || (~isgood(obs(i).tp(j).epKnh4))
                obs(i).tp(j).epKnh4        = nan;
                wobs(i,sys.tp(j).ipKnh4)   = (epKnh4).^(-2);
            else
                wobs(i,sys.tp(j).ipKnh4)   = (obs(i).tp(j).epKnh4).^(-2);
            end
            if (~isfield(obs(i).tp(j),'nh4')) || (~isgood(obs(i).tp(j).nh4))
                obs(i).tp(j).nh4           = nan;
                yobs(i,sys.tp(j).ipnh4)    = nan;
            else
                yobs(i,sys.tp(j).ipnh4)    = p(obs(i).tp(j).nh4*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'enh4')) || (~isgood(obs(i).tp(j).enh4))
                obs(i).tp(j).enh4          = nan;
                wobs(i,sys.tp(j).ipnh4)    = nan;
            else
                wobs(i,sys.tp(j).ipnh4) = w(obs(i).tp(j).nh4,obs(i).tp(j).enh4);
            end
            if (~isfield(obs(i).tp(j),'nh3')) || (~isgood(obs(i).tp(j).nh3))
                obs(i).tp(j).nh3           = nan;
                yobs(i,sys.tp(j).ipnh3)    = nan;
            else
                yobs(i,sys.tp(j).ipnh3)    = p(obs(i).tp(j).nh3*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'enh3')) || (~isgood(obs(i).tp(j).enh3))
                obs(i).tp(j).enh3          = nan;
                wobs(i,sys.tp(j).ipnh3)    = nan;
            else
                wobs(i,sys.tp(j).ipnh3) = w(obs(i).tp(j).nh3,obs(i).tp(j).enh3);
            end

            % sulfide system
            if (~isfield(obs(i).tp(j),'pKh2s')) || (~isgood(obs(i).tp(j).pKh2s))
                obs(i).tp(j).pKh2s         = nan;
                yobs(i,sys.tp(j).ipKh2s)   = pKh2s;
            else
                yobs(i,sys.tp(j).ipKh2s)   = obs(i).tp(j).pKh2s;
            end
            if (~isfield(obs(i).tp(j),'epKh2s')) || (~isgood(obs(i).tp(j).epKh2s))
                obs(i).tp(j).epKh2s        = nan;
                wobs(i,sys.tp(j).ipKh2s)   = (epKh2s).^(-2);
            else
                wobs(i,sys.tp(j).ipKh2s)   = (obs(i).tp(j).epKh2s).^(-2);
            end
            if (~isfield(obs(i).tp(j),'H2S')) || (~isgood(obs(i).tp(j).H2S))
                obs(i).tp(j).H2S           = nan;
                yobs(i,sys.tp(j).ipH2S)    = nan;
            else
                yobs(i,sys.tp(j).ipH2S) = p(obs(i).tp(j).H2S*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'eH2S')) || (~isgood(obs(i).tp(j).eH2S))
                obs(i).tp(j).eH2S          = nan;
                wobs(i,sys.tp(j).ipH2S)    = nan;
            else
                wobs(i,sys.tp(j).ipH2S) = w(obs(i).tp(j).H2S,obs(i).tp(j).eH2S);
            end            
            if (~isfield(obs(i).tp(j),'HS')) || (~isgood(obs(i).tp(j).HS))
                obs(i).tp(j).HS            = nan;
                yobs(i,sys.tp(j).ipHS)     = nan;
            else
                yobs(i,sys.tp(j).ipHS)     = p(obs(i).tp(j).HS*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'eHS')) || (~isgood(obs(i).tp(j).eHS))
                obs(i).tp(j).eHS           = nan;
                wobs(i,sys.tp(j).ipHS)     = nan;
            else
                wobs(i,sys.tp(j).ipHS) = w(obs(i).tp(j).HS,obs(i).tp(j).eHS);
            end
            
            % calcium carbonate solubility system
            if (~isfield(obs(i).tp(j),'pKar')) || (~isgood(obs(i).tp(j).pKar))
                obs(i).tp(j).pKar          = nan;
                yobs(i,sys.tp(j).ipKar)    = pKar;
            else
                yobs(i,sys.tp(j).ipKar)    = obs(i).tp(j).pKar;
            end
            if (~isfield(obs(i).tp(j),'epKar')) || (~isgood(obs(i).tp(j).epKar))
                obs(i).tp(j).epKar         = nan;
                wobs(i,sys.tp(j).ipKar)    = (epKar).^(-2);
            else
                wobs(i,sys.tp(j).ipKar)    = (obs(i).tp(j).epKar).^(-2);
            end
            if (~isfield(obs(i).tp(j),'pKca')) || (~isgood(obs(i).tp(j).pKca))
                obs(i).tp(j).pKca          = nan;
                yobs(i,sys.tp(j).ipKca)    = pKca;
            else
                yobs(i,sys.tp(j).ipKca)    = obs(i).tp(j).pKca;
            end
            if (~isfield(obs(i).tp(j),'epKca')) || (~isgood(obs(i).tp(j).epKca))
                obs(i).tp(j).epKca         = nan;
                wobs(i,sys.tp(j).ipKca)    = (epKca).^(-2);
            else
                wobs(i,sys.tp(j).ipKca)    = (obs(i).tp(j).epKca).^(-2);
            end
            if (~isfield(obs(i).tp(j),'OmegaAr')) || (~isgood(obs(i).tp(j).OmegaAr))
                obs(i).tp(j).OmegaAr           = nan;
                yobs(i,sys.tp(j).ipOmegaAr)    = nan;
            else
                yobs(i,sys.tp(j).ipOmegaAr) = p(obs(i).tp(j).OmegaAr); % Omega is dimensionless
            end
            if (~isfield(obs(i).tp(j),'eOmegaAr')) || (~isgood(obs(i).tp(j).eOmegaAr))
                obs(i).tp(j).eOmegaAr          = nan;
                wobs(i,sys.tp(j).ipOmegaAr)    = nan;
            else
                wobs(i,sys.tp(j).ipOmegaAr) = w(obs(i).tp(j).OmegaAr,obs(i).tp(j).eOmegaAr);
            end
            if (~isfield(obs(i).tp(j),'OmegaCa')) || (~isgood(obs(i).tp(j).OmegaCa))
                obs(i).tp(j).OmegaCa           = nan;
                yobs(i,sys.tp(j).ipOmegaCa)    = nan;
            else
                yobs(i,sys.tp(j).ipOmegaCa) = p(obs(i).tp(j).OmegaCa); % Omega is dimensionless
            end
            if (~isfield(obs(i).tp(j),'eOmegaCa')) || (~isgood(obs(i).tp(j).eOmegaCa))
                obs(i).tp(j).eOmegaCa          = nan;
                wobs(i,sys.tp(j).ipOmegaCa)    = nan;
            else
                wobs(i,sys.tp(j).ipOmegaCa) = w(obs(i).tp(j).OmegaCa,obs(i).tp(j).eOmegaCa);
            end
            if (~isfield(obs(i).tp(j),'ca')) || (~isgood(obs(i).tp(j).ca))
               obs(i).tp(j).ca        = nan;
                yobs(i,sys.tp(j).ipca) = nan; 
            else
                yobs(i,sys.tp(j).ipca) = p(obs(i).tp(j).ca*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'eca')) || (~isgood(obs(i).tp(j).eca))
                obs(i).tp(j).eca       = nan;
                wobs(i,sys.tp(j).ipca) = nan;
            else
                wobs(i,sys.tp(j).ipca)  = w(obs(i).tp(j).ca,obs(i).tp(j).eca);
            end

            if opt.pKalpha == 1
                % pKalpha system
                if (~isfield(obs(i).tp(j),'pKalpha')) || (~isgood(obs(i).tp(j).pKalpha))
                    obs(i).tp(j).pKalpha        = nan;
                    pKalpha = 4.0; % default
                    yobs(i,sys.tp(j).ipKalpha)  = pKalpha; 
                else
                    pKalpha = obs(i).tp(j).pKalpha;
                    yobs(i,sys.tp(j).ipKalpha)  = obs(i).tp(j).pKalpha;
                end
                if (~isfield(obs(i).tp(j),'epKalpha')) || (~isgood(obs(i).tp(j).epKalpha))
                    pKalpha = 4.0; % default
                    Kalpha     = q(pKalpha); % default
                    % pKalpha    = (4.0); % default
                    obs(i).tp(j).epKalpha       = nan;
                    wobs(i,sys.tp(j).ipKalpha)  = w(Kalpha,0.10*Kalpha); % 10%
                    % wobs(i,sys.tp(j).ipKalpha)  = (0.1*pKalpha)^(-2);
                else
                    wobs(i,sys.tp(j).ipKalpha)  = (obs(i).tp(j).pKalpha)^(-2);
                end
                if (pKalpha > 4.5) % > from Kerr, < from Humphreys
                    sys.M(sys.tp(j).row_alk, sys.tp(j).ipalpha) = 0;
                else
                    sys.M(sys.tp(j).row_alk, sys.tp(j).iphalpha) = 0;
                end
                % other unknowns in system
                obs(i).tp(j).palpha         = nan; % p(alpha)
                yobs(i,sys.tp(j).ipalpha)   = nan;
                obs(i).tp(j).epalpha        = nan;
                wobs(i,sys.tp(j).ipalpha)   = nan;
                obs(i).tp(j).phalpha        = nan; % p(H-alpha)
                yobs(i,sys.tp(j).iphalpha)  = nan;
                obs(i).tp(j).ephalpha       = nan;
                wobs(i,sys.tp(j).iphalpha)  = nan;
            end

            if opt.pKbeta == 1
                % pKbeta system
                if (~isfield(obs(i).tp(j),'pKbeta')) || (~isgood(obs(i).tp(j).pKbeta))
                    pKbeta                  = 6.95; % default
                    obs(i).tp(j).pKbeta         = nan;
                    yobs(i,sys.tp(j).ipKbeta)   = pKbeta; 
                else
                    pKbeta                      = obs(i).tp(j).pKbeta;
                    yobs(i,sys.tp(j).ipKbeta)   = obs(i).tp(j).pKbeta;
                end
                if (~isfield(obs(i).tp(j),'epKbeta')) || (~isgood(obs(i).tp(j).epKbeta))
                    pKbeta                      = 6.95; % default
                    Kbeta                       = q(pKbeta);
                    obs(i).tp(j).epKbeta        = nan;
                    wobs(i,sys.tp(j).ipKbeta)   = w(Kbeta,0.1*Kbeta);
                else
                    wobs(i,sys.tp(j).ipKbeta)   = (obs(i).tp(j).pKbeta)^(-2);
                end
                if (pKbeta > 4.5)
                    sys.M(sys.tp(j).row_alk, sys.tp(j).ipbeta) = 0;
                else
                    sys.M(sys.tp(j).row_alk, sys.tp(j).iphbeta) = 0;
                end
                % other unknowns in system
                obs(i).tp(j).pbeta          = nan; % p(beta)
                yobs(i,sys.tp(j).ipbeta)    = nan;
                obs(i).tp(j).epbeta         = nan;
                wobs(i,sys.tp(j).ipbeta)    = nan;
                obs(i).tp(j).phbeta         = nan; % p(H-beta)
                yobs(i,sys.tp(j).iphbeta)   = nan;
                obs(i).tp(j).ephbeta        = nan;
                wobs(i,sys.tp(j).iphbeta)   = nan;
            end

        end % for j = 1:nTP

    end
end



