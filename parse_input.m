function [obs,yobs,wobs] = parse_input(obs,sys,opt,nD)

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
    
        for j = 1:nTP % loop over (T,P) pairs
            
            % create obs structure with fieldnames
            if (~isfield(obs(i).tp(j),'epK0'))
                obs(i).tp(j).epK0   = [];
            end
            if (~isfield(obs(i).tp(j),'pK0'))
                obs(i).tp(j).pK0    = nan;
            end
            if (~isfield(obs(i).tp(j),'epK1'))
                obs(i).tp(j).epK1   = [];
            end
            if (~isfield(obs(i).tp(j),'pK1'))
                obs(i).tp(j).pK1    = nan;
            end
            if (~isfield(obs(i).tp(j),'epK2'))
                obs(i).tp(j).epK2   = [];
            end
            if (~isfield(obs(i).tp(j),'pK2'))
                obs(i).tp(j).pK2    = nan;
            end
            if (~isfield(obs(i).tp(j),'efco2')) || (~isgood(obs(i).tp(j).efco2))
                obs(i).tp(j).efco2  = [];
            end
            if (~isfield(obs(i).tp(j),'fco2')) || (~isgood(obs(i).tp(j).fco2))
                obs(i).tp(j).fco2   = nan;
            end
            if (~isfield(obs(i).tp(j),'eco2st')) || (~isgood(obs(i).tp(j).eco2st))
                obs(i).tp(j).eco2st = [];
            end
            if (~isfield(obs(i).tp(j),'co2st')) || (~isgood(obs(i).tp(j).co2st))
                obs(i).tp(j).co2st  = nan;
            end
            if (~isfield(obs(i).tp(j),'ehco3')) || (~isgood(obs(i).tp(j).ehco3))
                obs(i).tp(j).ehco3  = [];
            end
            if (~isfield(obs(i).tp(j),'hco3')) || (~isgood(obs(i).tp(j).hco3))
                obs(i).tp(j).hco3   = nan;
            end
            if (~isfield(obs(i).tp(j),'eco3')) || (~isgood(obs(i).tp(j).eco3))
                obs(i).tp(j).eco3   = [];
            end
            if (~isfield(obs(i).tp(j),'co3')) || (~isgood(obs(i).tp(j).co3))
                obs(i).tp(j).co3    = nan;
            end
            if (~isfield(obs(i).tp(j),'eph')) || (~isgood(obs(i).tp(j).eph))
                obs(i).tp(j).eph    = [];
            end
            if (~isfield(obs(i).tp(j),'ph')) || (~isgood(obs(i).tp(j).ph))
                obs(i).tp(j).ph     = nan;
            end

            % borate system
            if (~isfield(obs(i).tp(j),'epKb'))
                obs(i).tp(j).epKb   = [];
            end
            if (~isfield(obs(i).tp(j),'pKb'))
                obs(i).tp(j).pKb    = nan;
            end
            if (~isfield(obs(i).tp(j),'eboh4')) || (~isgood(obs(i).tp(j).eboh4))
                obs(i).tp(j).eboh4  = [];
            end
            if (~isfield(obs(i).tp(j),'boh4')) || (~isgood(obs(i).tp(j).boh4))
                obs(i).tp(j).boh4   = [];
            end
            if (~isfield(obs(i).tp(j),'eboh3')) || (~isgood(obs(i).tp(j).eboh3))
                obs(i).tp(j).eboh3  = [];
            end
            if (~isfield(obs(i).tp(j),'boh3')) || (~isgood(obs(i).tp(j).boh3))
                obs(i).tp(j).boh3   = [];
            end

            % water 
            if (~isfield(obs(i).tp(j), 'epKw'))
                obs(i).tp(j).epKw   = [];
            end
            if (~isfield(obs(i).tp(j), 'pKw'))
                obs(i).tp(j).pKw    = nan;
            end
            if (~isfield(obs(i).tp(j),'eoh')) || (~isgood(obs(i).tp(j).eoh))
                obs(i).tp(j).eoh    = [];
            end
            if (~isfield(obs(i).tp(j),'oh')) || (~isgood(obs(i).tp(j).oh))
                obs(i).tp(j).oh     = [];
            end
            
            % sulfate system
            if (~isfield(obs(i).tp(j), 'epKs'))
                obs(i).tp(j).epKs   = [];
            end
            if (~isfield(obs(i).tp(j), 'pKs'))
                obs(i).tp(j).pKs    = nan;
            end
            if (~isfield(obs(i).tp(j), 'so4')) || (~isgood(obs(i).tp(j).so4))
                obs(i).tp(j).so4    = [];
            end
            if (~isfield(obs(i).tp(j), 'eso4')) || (~isgood(obs(i).tp(j).eso4))
                obs(i).tp(j).eso4   = [];
            end
            if (~isfield(obs(i).tp(j), 'hso4')) || (~isgood(obs(i).tp(j).hso4))
                obs(i).tp(j).hso4   = [];
            end
            if (~isfield(obs(i).tp(j), 'ehso4')) || (~isgood(obs(i).tp(j).ehso4))
                obs(i).tp(j).ehso4  = [];
            end
            % fluoride system
            if (~isfield(obs(i).tp(j), 'epKf'))
                obs(i).tp(j).epKf   = [];
            end
            if (~isfield(obs(i).tp(j), 'pKf'))
                obs(i).tp(j).pKf    = nan;
            end
            if (~isfield(obs(i).tp(j), 'F')) || (~isgood(obs(i).tp(j).F))
                obs(i).tp(j).F      = [];
            end
            if (~isfield(obs(i).tp(j), 'eF')) || (~isgood(obs(i).tp(j).eF))
                obs(i).tp(j).eF     = [];
            end
            if (~isfield(obs(i).tp(j), 'HF')) || (~isgood(obs(i).tp(j).HF))
                obs(i).tp(j).HF     = [];
            end
            if (~isfield(obs(i).tp(j), 'eHF')) || (~isgood(obs(i).tp(j).eHF))
                obs(i).tp(j).eHF    = [];
            end
            % phosphate system
            if (~isfield(obs(i).tp(j), 'epKp1'))
                obs(i).tp(j).epKp1  = [];
            end
            if (~isfield(obs(i).tp(j), 'pKp1'))
                obs(i).tp(j).pKp1   = nan;
            end
            if (~isfield(obs(i).tp(j), 'epKp2'))
                obs(i).tp(j).epKp2  = [];
            end
            if (~isfield(obs(i).tp(j), 'pKp2'))
                obs(i).tp(j).pKp2   = nan;
            end
            if (~isfield(obs(i).tp(j), 'epKp3'))
                obs(i).tp(j).epKp3  = [];
            end
            if (~isfield(obs(i).tp(j), 'pKp3'))
                obs(i).tp(j).pKp3   = nan;
            end
            if (~isfield(obs(i).tp(j), 'h3po4')) || (~isgood(obs(i).tp(j).h3po4))
                obs(i).tp(j).h3po4  = [];
            end
            if (~isfield(obs(i).tp(j), 'eh3po4')) || (~isgood(obs(i).tp(j).eh3po4))
                obs(i).tp(j).eh3po4 = [];
            end
            if (~isfield(obs(i).tp(j), 'h2po4')) || (~isgood(obs(i).tp(j).h2po4))
                obs(i).tp(j).h2po4  = [];
            end
            if (~isfield(obs(i).tp(j), 'eh2po4')) || (~isgood(obs(i).tp(j).eh2po4))
                obs(i).tp(j).eh2po4 = [];
            end
            if (~isfield(obs(i).tp(j), 'hpo4')) || (~isgood(obs(i).tp(j).hpo4))
                obs(i).tp(j).hpo4   = [];
            end
            if (~isfield(obs(i).tp(j), 'ehpo4')) || (~isgood(obs(i).tp(j).ehpo4))
                obs(i).tp(j).ehpo4  = [];
            end
            if (~isfield(obs(i).tp(j), 'po4')) || (~isgood(obs(i).tp(j).po4))
                obs(i).tp(j).po4    = [];
            end
            if (~isfield(obs(i).tp(j), 'epo4')) || (~isgood(obs(i).tp(j).epo4))
                obs(i).tp(j).epo4   = [];
            end
            
            if (~isfield(obs(i).tp(j), 'epKsi'))
                obs(i).tp(j).epKsi  = [];
            end
            if (~isfield(obs(i).tp(j), 'pKsi'))
                obs(i).tp(j).pKsi   = nan;
            end
            if (~isfield(obs(i).tp(j), 'sioh4')) || (~isgood(obs(i).tp(j).sioh4))
                obs(i).tp(j).sioh4      = [];
            end
            if (~isfield(obs(i).tp(j), 'esioh4')) || (~isgood(obs(i).tp(j).esioh4))
                obs(i).tp(j).esioh4     = [];
            end
            if (~isfield(obs(i).tp(j), 'siooh3')) || (~isgood(obs(i).tp(j).siooh3))
                obs(i).tp(j).siooh3     = [];
            end
            if (~isfield(obs(i).tp(j), 'esiooh3')) || (~isgood(obs(i).tp(j).esiooh3))
                obs(i).tp(j).esiooh3    = [];
            end
            if (~isfield(obs(i).tp(j), 'epKnh4'))
                obs(i).tp(j).epKnh4     = [];
            end
            if (~isfield(obs(i).tp(j), 'pKnh4'))
                obs(i).tp(j).pKnh4      = nan;
            end
            if (~isfield(obs(i).tp(j), 'nh3')) || (~isgood(obs(i).tp(j).nh3))
                obs(i).tp(j).nh3    = [];
            end
            if (~isfield(obs(i).tp(j), 'enh3')) || (~isgood(obs(i).tp(j).enh3))
                obs(i).tp(j).enh3   = [];
            end
            if (~isfield(obs(i).tp(j), 'nh4')) || (~isgood(obs(i).tp(j).nh4))
                obs(i).tp(j).nh4    = [];
            end
            if (~isfield(obs(i).tp(j), 'enh4')) || (~isgood(obs(i).tp(j).enh4))
                obs(i).tp(j).enh4   = [];
            end
            if (~isfield(obs(i).tp(j), 'epKh2s'))
                obs(i).tp(j).epKh2s = [];
            end
            if (~isfield(obs(i).tp(j), 'pKh2s'))
                obs(i).tp(j).pKh2s  = nan;
            end
            if (~isfield(obs(i).tp(j), 'HS')) || (~isgood(obs(i).tp(j).HS))
                obs(i).tp(j).HS     = [];
            end
            if (~isfield(obs(i).tp(j), 'eHS')) || (~isgood(obs(i).tp(j).eHS))
                obs(i).tp(j).eHS    = [];
            end
            if (~isfield(obs(i).tp(j), 'H2S')) || (~isgood(obs(i).tp(j).H2S))
                obs(i).tp(j).H2S    = [];
            end
            if (~isfield(obs(i).tp(j), 'eH2S')) || (~isgood(obs(i).tp(j).eH2S))
                obs(i).tp(j).eH2S   = [];
            end
            if (~isfield(obs(i).tp(j),'epco2')) || (~isgood(obs(i).tp(j).epco2))
                obs(i).tp(j).epco2  = [];
            end
            if (~isfield(obs(i).tp(j),'epp2f')) || (~isgood(obs(i).tp(j).epp2f))
                obs(i).tp(j).epp2f  = [];
            end
            if (~isfield(obs(i).tp(j),'pp2f')) || (~isgood(obs(i).tp(j).pp2f))
                obs(i).tp(j).pp2f   = nan;
            end
            if (~isfield(obs(i).tp(j),'pco2')) || (~isgood(obs(i).tp(j).pco2))
                obs(i).tp(j).pco2   = nan;
            end

            if (~isfield(obs(i).tp(j), 'epKar'))
                obs(i).tp(j).epKar = [];
            end
            if (~isfield(obs(i).tp(j), 'pKar'))
                obs(i).tp(j).pKar   = nan;
            end
            if (~isfield(obs(i).tp(j), 'ca')) || (~isgood(obs(i).tp(j).ca))
                obs(i).tp(j).ca     = [];
            end
            if (~isfield(obs(i).tp(j), 'eca')) || (~isgood(obs(i).tp(j).eca))
                obs(i).tp(j).eca    = [];
            end
            if (~isfield(obs(i).tp(j), 'OmegaAr')) || (~isgood(obs(i).tp(j).OmegaAr))
                obs(i).tp(j).OmegaAr    = [];
            end
            if (~isfield(obs(i).tp(j), 'eOmegaAr')) || (~isgood(obs(i).tp(j).eOmegaAr))
                obs(i).tp(j).eOmegaAr   = [];
            end 
            if (~isfield(obs(i).tp(j), 'pKar'))
                obs(i).tp(j).pKar       = nan;
            end
            if (~isfield(obs(i).tp(j), 'epKca'))
                obs(i).tp(j).epKar      = [];
            end
            if (~isfield(obs(i).tp(j), 'pKca'))
                obs(i).tp(j).pKca       = nan;
            end
            if (~isfield(obs(i).tp(j), 'epKca'))
                obs(i).tp(j).epKca      = [];
            end
            if (~isfield(obs(i).tp(j), 'OmegaCa')) || (~isgood(obs(i).tp(j).OmegaCa))
                obs(i).tp(j).OmegaCa    = [];
            end
            if (~isfield(obs(i).tp(j), 'eOmegaCa')) || (~isgood(obs(i).tp(j).eOmegaCa))
                obs(i).tp(j).eOmegaCa   = [];
            end
            
            if (~isfield(obs(i).tp(j),'eph_free')) || (~isgood(obs(i).tp(j).eph_free))
                obs(i).tp(j).eph_free = [];
            end
            if (~isfield(obs(i).tp(j),'ph_free')) || (~isgood(obs(i).tp(j).ph_free))
                obs(i).tp(j).ph_free = nan;
            end
            if (~isfield(obs(i).tp(j),'epfH')) || (~isgood(obs(i).tp(j).epfH))
                obs(i).tp(j).epfH   = [];
            end
            if (~isfield(obs(i).tp(j),'pfH')) || (~isgood(obs(i).tp(j).pfH))
                obs(i).tp(j).pfH    = nan;
            end
        end

        for ii = 1:nTP % loop over all the pressure and temperature sensitive components
            yobs(i,sys.tp(ii).iT)   = obs(i).tp(ii).T;
            yobs(i,sys.tp(ii).iP)   = obs(i).tp(ii).P;
            wobs(i,sys.tp(ii).iT)   = (obs(i).tp(ii).eT)^(-2);
            wobs(i,sys.tp(ii).iP)   = (obs(i).tp(ii).eP)^(-2);
            [pK,~,epK]              = calc_pK(opt, obs(i).tp(ii).T, ...      % T,
                                              obs(i).sal, obs(i).tp(ii).P ); % S, P
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
            
            %
            % co2 solubility and fugacity
            %
            if (isgood(obs(i).tp(ii).epK0))
                wobs(i,sys.tp(ii).ipK0) = (obs(i).tp(ii).epK0)^(-2);
            else
                obs(i).tp(ii).epK0 = epK0;
                wobs(i,sys.tp(ii).ipK0) = (obs(i).tp(ii).epK0)^(-2);
            end
            if (isgood(obs(i).tp(ii).pK0))
                yobs(i,sys.tp(ii).ipK0) = obs(i).tp(ii).pK0;
            else
                yobs(i,sys.tp(ii).ipK0) = pK0;
                obs(i).tp(ii).pK0 = nan;
            end
            if (isgood(obs(i).tp(ii).fco2))
                yobs(i,sys.tp(ii).ipfco2) = p(obs(i).tp(ii).fco2*1e-6); % convt µatm to atm
            else
                yobs(i,sys.tp(ii).ipfco2) = nan;
                obs(i).tp(ii).fco2 = nan;
            end
            if (isgood(obs(i).tp(ii).efco2))
                wobs(i,sys.tp(ii).ifco2) = w(obs(i).tp(ii).fco2, obs(i).tp(ii).efco2);
            else
                wobs(i,sys.tp(ii).ipfco2) = nan;
                obs(i).tp(ii).efco2 = nan;
            end
            %
            % carbonate system
            %
            if (isgood(obs(i).tp(ii).epK1))
                wobs(i,sys.tp(ii).ipK1) = (obs(i).tp(ii).epK1)^(-2);
            else
                obs(i).tp(ii).epK1 = epK1;
                wobs(i,sys.tp(ii).ipK1) = (obs(i).tp(ii).epK1)^(-2);
            end
            if (isgood(obs(i).tp(ii).pK1))
                yobs(i,sys.tp(ii).ipK1) = obs(i).tp(ii).pK1;
            else
                yobs(i,sys.tp(ii).ipK1) = pK1;
                obs(i).tp(ii).pK1  = nan;
            end
            if (isgood(obs(i).tp(ii).epK2))
                wobs(i,sys.tp(ii).ipK2) = (obs(i).tp(ii).epK2)^(-2);
            else
                obs(i).tp(ii).epK2 = epK2;
                wobs(i,sys.tp(ii).ipK2)  = (obs(i).tp(ii).epK2)^(-2);
            end
            if (isgood(obs(i).tp(ii).pK2))
                yobs(i,sys.tp(ii).ipK2) = obs(i).tp(ii).pK2;
            else
                yobs(i,sys.tp(ii).ipK2) = pK2;
                obs(i).tp(ii).pK2 = nan;
            end
            if (isgood(obs(i).tp(ii).co2st))
                yobs(i,sys.tp(ii).ipco2st) = p(obs(i).tp(ii).co2st*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipco2st) = nan;
                obs(i).tp(ii).co2st = nan;
            end
            if (isgood(obs(i).tp(ii).eco2st))
                wobs(i,sys.tp(ii).ipco2st) = w(obs(i).tp(ii).co2st, ...
                                               obs(i).tp(ii).eco2st);
            else
                wobs(i,sys.tp(ii).ipco2st) = nan;
                obs(i).tp(ii).eco2st = nan;
            end
            if (isgood(obs(i).tp(ii).hco3))
                yobs(i,sys.tp(ii).iphco3) = p(obs(i).tp(ii).hco3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).iphco3) = nan;
                obs(i).tp(ii).hco3 = nan;
            end
            if (isgood(obs(i).tp(ii).ehco3))
                wobs(i,sys.tp(ii).iphco3) = w(obs(i).tp(ii).hco3, obs(i).tp(ii).ehco3);
            else
                wobs(i,sys.tp(ii).iphco3) = nan;
                obs(i).tp(ii).ehco3 = nan;
            end
            if (isgood(obs(i).tp(ii).co3))
                yobs(i,sys.tp(ii).ipco3) = p(obs(i).tp(ii).co3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipco3) = nan;
                obs(i).tp(ii).co3 = nan;
            end
            if (isgood(obs(i).tp(ii).eco3))
                wobs(i,sys.tp(ii).ipco3) = w(obs(i).tp(ii).co3, obs(i).tp(ii).eco3);
            else
                wobs(i,sys.tp(ii).ipco3) = nan;
                obs(i).tp(ii).eco3 = nan;
            end
            if (isgood(obs(i).tp(ii).ph))
                yobs(i,sys.tp(ii).iph) = obs(i).tp(ii).ph ;
            else
                yobs(i,sys.tp(ii).iph) = nan;
                obs(i).tp(ii).ph = nan;
            end
            if (isgood(obs(i).tp(ii).eph))
                wobs(i,sys.tp(ii).iph) = (obs(i).tp(ii).eph).^(-2);
            else
                wobs(i,sys.tp(ii).iph) = nan;
                obs(i).tp(ii).eph = nan;
            end
            % borate system 
            if (isgood(obs(i).tp(ii).epKb))
                wobs(i,sys.tp(ii).ipKb) = (obs(i).tp(ii).epKb).^(-2);
            else
                obs(i).tp(ii).epKb = epKb;
                wobs(i,sys.tp(ii).ipKb)  = (obs(i).tp(ii).epKb).^(-2);  % wKb = 1/(1 + (0.01/pKsys(4)))^2 ;
            end
            if (isgood(obs(i).tp(ii).pKb))
                yobs(i,sys.tp(ii).ipKb) = obs(i).tp(ii).pKb;
            else
                yobs(i,sys.tp(ii).ipKb) = pKb; 
                obs(i).tp(ii).pKb = nan;
            end
            if (isgood(obs(i).tp(ii).boh3))
                yobs(i,sys.tp(ii).ipboh3) = p(obs(i).tp(ii).boh3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipboh3) = nan;
                obs(i).tp(ii).boh3 = nan;
            end
            if (isgood(obs(i).tp(ii).eboh3))
                wobs(i,sys.tp(ii).ipboh3) = w(obs(i).tp(ii).boh3, obs(i).tp(ii).eboh3);
            else
                wobs(i,sys.tp(ii).ipboh3) = nan;
                obs(i).tp(ii).eboh3 = nan;
            end
            if (isgood(obs(i).tp(ii).boh4))
                yobs(i,sys.tp(ii).ipboh4) = p(obs(i).tp(ii).boh4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipboh4) = nan;
                obs(i).tp(ii).boh4 = nan;
            end
            if (isgood(obs(i).tp(ii).eboh4))
                wobs(i,sys.tp(ii).ipboh4) = w(obs(i).tp(ii).boh4, obs(i).tp(ii).eboh4);
            else
                wobs(i,sys.tp(ii).ipboh4) = nan;
                obs(i).tp(ii).eboh4 = nan;
            end
            
            % water dissociation
            if (isgood(obs(i).tp(ii).epKw))
                wobs(i,sys.tp(ii).ipKw) = (obs(i).tp(ii).epKw).^(-2);
            else
                obs(i).tp(ii).epKw = epKw;
                wobs(i,sys.tp(ii).ipKw) = (obs(i).tp(ii).epKw).^(-2);
            end
            if (isgood(obs(i).tp(ii).pKw))
                yobs(i,sys.tp(ii).ipKw) = obs(i).tp(ii).pKw;
            else
                yobs(i,sys.tp(ii).ipKw) = pKw;
                obs(i).tp(ii).pKw = nan;
            end
            if (isgood(obs(i).tp(ii).oh))
                yobs(i,sys.tp(ii).ipoh) = p(obs(i).tp(ii).oh*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipoh) = nan;
                obs(i).tp(ii).oh = nan;
            end
            if (isgood(obs(i).tp(ii).eoh))
                wobs(i,sys.tp(ii).ipoh) = w(obs(i).tp(ii).oh, obs(i).tp(ii).eoh);
            else
                wobs(i,sys.tp(ii).ipoh) = nan;
                obs(i).tp(ii).eoh = nan;
            end
            
            % sulfate system
            if (isgood(obs(i).tp(ii).epKs))
                wobs(i,sys.tp(ii).ipKs) = (obs(i).tp(ii).epKs).^(-2);
            else
                obs(i).tp(ii).epKs = epKs;
                wobs(i,sys.tp(ii).ipKs) = (obs(i).tp(ii).epKs).^(-2); % wKs = 1/(1 + (0.0021/pKsys(6)))^2 ;
            end
            if (isgood(obs(i).tp(ii).pKs))
                yobs(i,sys.tp(ii).ipKs) = obs(i).tp(ii).pKs;
            else
                yobs(i,sys.tp(ii).ipKs) = pKs;
                obs(i).tp(ii).pKs = nan;
            end
            if (isgood(obs(i).tp(ii).hso4))
                yobs(i,sys.tp(ii).iphso4) = p(obs(i).tp(ii).hso4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).iphso4) = nan;
                obs(i).tp(ii).hso4 = nan;
            end
            if (isgood(obs(i).tp(ii).ehso4))
                wobs(i,sys.tp(ii).iphso4) = w(obs(i).tp(ii).hso4, obs(i).tp(ii).ehso4);
            else
                wobs(i,sys.tp(ii).iphso4) = nan;
                obs(i).tp(ii).ehso4 = nan;
            end
            if (isgood(obs(i).tp(ii).so4))
                yobs(i,sys.tp(ii).ipso4) = p(obs(i).tp(ii).so4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipso4) = nan;
                obs(i).tp(ii).so4 = nan;
            end
            if (isgood(obs(i).tp(ii).eso4))
                wobs(i,sys.tp(ii).ipso4) = w(obs(i).tp(ii).so4, obs(i).tp(ii).eso4);
            else
                wobs(i,sys.tp(ii).ipso4) = nan;
                obs(i).tp(ii).eso4 = nan;
            end
            % fluoride system
            if (isgood(obs(i).tp(ii).epKf))
                wobs(i,sys.tp(ii).ipKf) = (obs(i).tp(ii).epKf).^(-2);
            else
                obs(i).tp(ii).epKf = epKf;
                wobs(i,sys.tp(ii).ipKf) = (obs(i).tp(ii).epKf).^(-2);
            end
            if (isgood(obs(i).tp(ii).pKf))
                yobs(i,sys.tp(ii).ipKf) = obs(i).tp(ii).pKf;
            else
                yobs(i,sys.tp(ii).ipKf) = pKf;
                obs(i).tp(ii).pKf = nan;
            end
            if (isgood(obs(i).tp(ii).F))
                yobs(i,sys.tp(ii).ipF) = p(obs(i).tp(ii).F*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipF) = nan;
                obs(i).tp(ii).F = nan;
            end
            if (isgood(obs(i).tp(ii).eF))
                wobs(i,sys.tp(ii).ipF) = w(obs(i).tp(ii).F, obs(i).tp(ii).eF);
            else
                wobs(i,sys.tp(ii).ipF) = nan;
                obs(i).tp(ii).eF = nan;
            end
            if (isgood(obs(i).tp(ii).HF))
                yobs(i,sys.tp(ii).ipHF) = p(obs(i).tp(ii).HF*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipHF) = nan;
                obs(i).tp(ii).HF = nan;
            end
            if (isgood(obs(i).tp(ii).eHF))
                wobs(i,sys.tp(ii).ipHF) = w(obs(i).tp(ii).HF, obs(i).tp(ii).eHF);
            else
                wobs(i,sys.tp(ii).ipHF) = nan;
                obs(i).tp(ii).eHF = nan;
            end
            
            % phosphate system
            if (isgood(obs(i).tp(ii).epKp1))
                wobs(i,sys.tp(ii).ipKp1) = (obs(i).tp(ii).epKp1).^(-2);
            else
                obs(i).tp(ii).epKp1 = epKp1;
                wobs(i,sys.tp(ii).ipKp1) = (obs(i).tp(ii).epKp1).^(-2);
            end
            if (isgood(obs(i).tp(ii).pKp1))
                yobs(i,sys.tp(ii).ipKp1) = obs(i).tp(ii).pKp1;
            else
                yobs(i,sys.tp(ii).ipKp1) = pKp1;
                obs(i).tp(ii).pKp1 = nan;
            end
            if (isgood(obs(i).tp(ii).epKp2))
                wobs(i,sys.tp(ii).ipKp2) = (obs(i).tp(ii).epKp2).^(-2);
            else
                obs(i).tp(ii).epKp2 = epKp2;
                wobs(i,sys.tp(ii).ipKp2) = (obs(i).tp(ii).epKp2).^(-2); 
            end
            if (isgood(obs(i).tp(ii).pKp2))
                yobs(i,sys.tp(ii).ipKp2) = obs(i).tp(ii).pKp2;
            else
                yobs(i,sys.tp(ii).ipKp2) = pKp2;
                obs(i).tp(ii).pKp2 = nan;
            end
            if (isgood(obs(i).tp(ii).epKp3))
                wobs(i,sys.tp(ii).ipKp3) = (obs(i).tp(ii).epKp3).^(-2);
            else
                obs(i).tp(ii).epKp3 = epKp3;
                wobs(i,sys.tp(ii).ipKp3) = (obs(i).tp(ii).epKp3).^(-2);
            end
            if (isgood(obs(i).tp(ii).pKp3))
                yobs(i,sys.tp(ii).ipKp3) = obs(i).tp(ii).pKp3;
            else
                yobs(i,sys.tp(ii).ipKp3) = pKp3;
                obs(i).tp(ii).pKp3 = nan;
            end
            if (isgood(obs(i).tp(ii).h3po4))
                yobs(i,sys.tp(ii).iph3po4) = p(obs(i).tp(ii).h3po4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).iph3po4) = nan;
                obs(i).tp(ii).h3po4 = nan;
            end
            if (isgood(obs(i).tp(ii).eh3po4))
                wobs(i,sys.tp(ii).iph3po4) = w(obs(i).tp(ii).h3po4, obs(i).tp(ii).eh3po4);
            else
                wobs(i,sys.tp(ii).iph3po4) = nan;
                obs(i).tp(ii).eh3po4 = nan;
            end
            if (isgood(obs(i).tp(ii).h2po4))
                yobs(i,sys.tp(ii).iph2po4) = p(obs(i).tp(ii).h2po4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).iph2po4) = nan;
                obs(i).tp(ii).h2po4 = nan;
            end
            if (isgood(obs(i).tp(ii).eh2po4))
                wobs(i,sys.tp(ii).iph2po4) = w(obs(i).tp(ii).h2po4, obs(i).tp(ii).eh2po4);
            else
                wobs(i,sys.tp(ii).iph2po4) = nan;
                obs(i).tp(ii).eh2po4 = nan;
            end
            if (isgood(obs(i).tp(ii).hpo4))
                yobs(i,sys.tp(ii).iphpo4) = p(obs(i).tp(ii).hpo4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).iphpo4) = nan;
                obs(i).tp(ii).hpo4 = nan;
            end
            if (isgood(obs(i).tp(ii).ehpo4))
                wobs(i,sys.tp(ii).iphpo4) = w(obs(i).tp(ii).hpo4, obs(i).tp(ii).ehpo4);
            else
                wobs(i,sys.tp(ii).iphpo4) = nan;
                obs(i).tp(ii).ehpo4 = nan;
            end
            if (isgood(obs(i).tp(ii).po4))
                yobs(i,sys.tp(ii).ippo4) = p(obs(i).tp(ii).po4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ippo4) = nan;
                obs(i).tp(ii).po4 = nan;
            end
            if (isgood(obs(i).tp(ii).epo4))
                wobs(i,sys.tp(ii).ippo4) = w(obs(i).tp(ii).po4, obs(i).tp(ii).epo4);
            else
                wobs(i,sys.tp(ii).ippo4) = nan;
                obs(i).tp(ii).epo4 = nan;
            end
            % silicate system
            if (isgood(obs(i).tp(ii).epKsi))
                wobs(i,sys.tp(ii).ipKsi) = (obs(i).tp(ii).epKsi).^(-2);
            else
                obs(i).tp(ii).epKsi = epKsi;
                wobs(i,sys.tp(ii).ipKsi) = (obs(i).tp(ii).epKsi).^(-2);  % wKSi = 1/(1 + (0.02/pKsys(11)))^2 ;
            end
            if (isgood(obs(i).tp(ii).pKsi))
                yobs(i,sys.tp(ii).ipKsi) = obs(i).tp(ii).pKsi;
            else
                yobs(i,sys.tp(ii).ipKsi) = pKsi;
                obs(i).tp(ii).pKsi = nan;
            end
            if (isgood(obs(i).tp(ii).sioh4))
                yobs(i,sys.tp(ii).ipsioh4) = p(obs(i).tp(ii).sioh4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipsioh4) = nan;
                obs(i).tp(ii).sioh4 = nan;
            end
            if (isgood(obs(i).tp(ii).esioh4))
                wobs(i,sys.tp(ii).ipsioh4) = w(obs(i).tp(ii).sioh4, obs(i).tp(ii).esioh4);
            else
                wobs(i,sys.tp(ii).ipsioh4) = nan;
                obs(i).tp(ii).esioh4 = nan;
            end
            if (isgood(obs(i).tp(ii).siooh3))
                yobs(i,sys.tp(ii).ipsiooh3) = p(obs(i).tp(ii).siooh3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipsiooh3) = nan;
                obs(i).tp(ii).siooh3 = nan;
            end
            if (isgood(obs(i).tp(ii).esiooh3))
                wobs(i,sys.tp(ii).ipsiooh3) = w(obs(i).tp(ii).siooh3, obs(i).tp(ii).esiooh3);
            else
                wobs(i,sys.tp(ii).ipsiooh3) = nan;
                obs(i).tp(ii).esiooh3 = nan;
            end
            % amonia system
            if (isgood(obs(i).tp(ii).epKnh4))
                wobs(i,sys.tp(ii).ipKnh4) = (obs(i).tp(ii).epKnh4).^(-2);
            else
                obs(i).tp(ii).epKnh4 = epKnh4;
                wobs(i,sys.tp(ii).ipKnh4) = (obs(i).tp(ii).epKnh4).^(-2);  % wKnh4 = 1/(1 + (0.00017/pKsys(11)))^2 ;
            end
            if (isgood(obs(i).tp(ii).pKnh4))
                yobs(i,sys.tp(ii).ipKnh4) = obs(i).tp(ii).pKnh4;
            else
                yobs(i,sys.tp(ii).ipKnh4) = pKnh4;
                obs(i).tp(ii).pKnh4 = nan;
            end
            if (isgood(obs(i).tp(ii).nh4))
                yobs(i,sys.tp(ii).ipnh4) = p(obs(i).tp(ii).nh4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipnh4) = nan;
                obs(i).tp(ii).nh4 = nan;
            end
            if (isgood(obs(i).tp(ii).enh4))
                wobs(i,sys.tp(ii).ipnh4) = w(obs(i).tp(ii).nh4, obs(i).tp(ii).enh4);
            else
                wobs(i,sys.tp(ii).ipnh4) = nan;
                obs(i).tp(ii).enh4 = nan;
            end
            if (isgood(obs(i).tp(ii).nh3))
                yobs(i,sys.tp(ii).ipnh3) = p(obs(i).tp(ii).nh3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipnh3) = nan;
                obs(i).tp(ii).nh3 = nan;
            end
            if (isgood(obs(i).tp(ii).enh3))
                wobs(i,sys.tp(ii).ipnh3) = w(obs(i).tp(ii).nh3, obs(i).tp(ii).enh3);
            else
                wobs(i,sys.tp(ii).ipnh3) = nan;
                obs(i).tp(ii).enh3 = nan;
            end
            % sulfide system
            if (isgood(obs(i).tp(ii).epKh2s))
                wobs(i,sys.tp(ii).ipKh2s) = (obs(i).tp(ii).epKh2s).^(-2);
            else
                obs(i).tp(ii).epKh2s = epKh2s;
                wobs(i,sys.tp(ii).ipKh2s) = (obs(i).tp(ii).epKh2s).^(-2);  % wKh2s = 1/(1 + (0.033/pKsys(11)))^2 ;
            end
            if (isgood(obs(i).tp(ii).pKh2s))
                yobs(i,sys.tp(ii).ipKh2s) = obs(i).tp(ii).pKh2s;
            else
                yobs(i,sys.tp(ii).ipKh2s) = pKh2s;
                obs(i).tp(ii).pKh2s = nan;
            end
            if (isgood(obs(i).tp(ii).H2S))
                yobs(i,sys.tp(ii).ipH2S) = p(obs(i).tp(ii).H2S*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipH2S) = nan;
                obs(i).tp(ii).H2S = nan;
            end
            if (isgood(obs(i).tp(ii).eH2S))
                wobs(i,sys.tp(ii).ipH2S) = w(obs(i).tp(ii).H2S, obs(i).tp(ii).eH2S);
            else
                wobs(i,sys.tp(ii).ipH2S) = nan;
                obs(i).tp(ii).eH2S = nan;
            end
            if (isgood(obs(i).tp(ii).HS))
                yobs(i,sys.tp(ii).ipHS) = p(obs(i).tp(ii).HS*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipHS) = nan;
                obs(i).tp(ii).HS = nan;
            end
            if (isgood(obs(i).tp(ii).eHS))
                wobs(i,sys.tp(ii).ipHS) = w(obs(i).tp(ii).HS, obs(i).tp(ii).eHS);
            else
                wobs(i,sys.tp(ii).ipHS) = nan;
                obs(i).tp(ii).eHS = nan;
            end
            if (isgood(obs(i).tp(ii).epp2f))
                wobs(i,sys.tp(ii).ipp2f) = (obs(i).tp(ii).epp2f)^(-2);
            else
                obs(i).tp(ii).epp2f = epp2f;
                wobs(i,sys.tp(ii).ipp2f) = (obs(i).tp(ii).epp2f)^(-2);
            end
            if (isgood(obs(i).tp(ii).pp2f))
                yobs(i,sys.tp(ii).ipp2f) = obs(i).tp(ii).pp2f;
            else
                yobs(i,sys.tp(ii).ipp2f) = pp2f;
                obs(i).tp(ii).pp2f = nan;
            end
            if (isgood(obs(i).tp(ii).pco2))
                yobs(i,sys.tp(ii).ippco2) = p(obs(i).tp(ii).pco2*1e-6); % convt µatm to atm
            else
                yobs(i,sys.tp(ii).ippco2) = nan;
                obs(i).tp(ii).pco2 = nan;
            end
            if (isgood(obs(i).tp(ii).epco2))
                wobs(i,sys.tp(ii).ippco2) = w(obs(i).tp(ii).pco2, obs(i).tp(ii).epco2);
            else
                wobs(i,sys.tp(ii).ippco2) = nan;
                obs(i).tp(ii).epco2 = nan;
            end
            
            % calcium carbonate mineral solubility
            if (isgood(obs(i).tp(ii).epKar))
                wobs(i,sys.tp(ii).ipKar) = (obs(i).tp(ii).epKar).^(-2);
            else
                obs(i).tp(ii).epKar = epKar;
                wobs(i,sys.tp(ii).ipKar) = (obs(i).tp(ii).epKar).^(-2);
            end
            if (isgood(obs(i).tp(ii).pKar))
                yobs(i,sys.tp(ii).ipKar) = obs(i).tp(ii).pKar;
            else
                yobs(i,sys.tp(ii).ipKar) = pKar;
                obs(i).tp(ii).pKar = nan;
            end
            if (isgood(obs(i).tp(ii).ca))
                yobs(i,sys.tp(ii).ipca) = p(obs(i).tp(ii).ca*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipca) = nan;
                obs(i).tp(ii).ca = nan;
            end
            if (isgood(obs(i).tp(ii).eca))
                wobs(i,sys.tp(ii).ipca)  = w(obs(i).tp(ii).ca, obs(i).tp(ii).eca);
            else
                wobs(i,sys.tp(ii).ipca) = nan;
                obs(i).tp(ii).eca = nan;
            end
            if (isgood(obs(i).tp(ii).OmegaAr))
                yobs(i,sys.tp(ii).ipOmegaAr) = p(obs(i).tp(ii).OmegaAr); % Omega is dimensionless
            else
                yobs(i,sys.tp(ii).ipOmegaAr) = nan;
                obs(i).tp(ii).OmegaAr = nan;
            end
            if (isgood(obs(i).tp(ii).eOmegaAr))
                wobs(i,sys.tp(ii).ipOmegaAr) = w(obs(i).tp(ii).OmegaAr, obs(i).tp(ii).eOmegaAr);
            else
                wobs(i,sys.tp(ii).ipOmegaAr) = nan;
                obs(i).tp(ii).eOmegaAr = nan;
            end
            if (isgood(obs(i).tp(ii).epKca))
                wobs(i,sys.tp(ii).ipKca) = (obs(i).tp(ii).epKca).^(-2);
            else
                obs(i).tp(ii).epKca = epKca;
                wobs(i,sys.tp(ii).ipKca) = (obs(i).tp(ii).epKca).^(-2);
            end
            if (isgood(obs(i).tp(ii).pKca))
                yobs(i,sys.tp(ii).ipKca) = obs(i).tp(ii).pKca;
            else
                yobs(i,sys.tp(ii).ipKca) = pKca;
                obs(i).tp(ii).pKca = nan;
            end
            if (isgood(obs(i).tp(ii).OmegaCa))
                yobs(i,sys.tp(ii).ipOmegaCa) = p(obs(i).tp(ii).OmegaCa); % Omega is dimensionless
            else
                yobs(i,sys.tp(ii).ipOmegaCa) = nan;
                obs(i).tp(ii).OmegaCa = nan;
            end
            if (isgood(obs(i).tp(ii).eOmegaCa))
                wobs(i,sys.tp(ii).ipOmegaCa) = w(obs(i).tp(ii).OmegaCa, obs(i).tp(ii).eOmegaCa);
            else
                wobs(i,sys.tp(ii).ipOmegaCa) = nan;
                obs(i).tp(ii).eOmegaCa = nan;
            end
            
            if (isgood(obs(i).tp(ii).ph_free)) % ph_free (ph on free scale)
                yobs(i,sys.tp(ii).iph_free) = obs(i).tp(ii).ph_free ;
            else
                yobs(i,sys.tp(ii).iph_free) = nan;
                obs(i).tp(ii).ph_free = nan;
            end
            if (isgood(obs(i).tp(ii).eph_free)) 
                wobs(i,sys.tp(ii).iph_free) = obs(i).tp(ii).eph_free.^(-2);
            else
                wobs(i,sys.tp(ii).iph_free) = nan;
                obs(i).tp(ii).eph_free = nan;
            end
            if (isgood(obs(i).tp(ii).pfH)) % pfH activity coefficient
                yobs(i,sys.tp(ii).ipfH) = obs(i).tp(ii).pfH ;
            else
                yobs(i,sys.tp(ii).ipfH) = pfH;
                obs(i).tp(ii).pfH = pfH;
            end
            if (isgood(obs(i).tp(ii).epfH)) % pfH activity coefficient
                wobs(i,sys.tp(ii).ipfH) = (obs(i).tp(ii).epfH).^(-2) ;
            else
                wobs(i,sys.tp(ii).ipfH) = (epfH).^(-2);
                obs(i).tp(ii).epfH = epfH;
            end
        end
    end
end



