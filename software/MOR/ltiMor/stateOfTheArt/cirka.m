function [sysr, s0, kIter, kIrkaTot, nSysm] = cirka(sys, s0, Opts) 

    %% Define execution options
    Def.qm0     = length(s0)+2;
    Def.s0m     = zeros(1,Def.qm0); 
    Def.maxiter = 8; Def.tol = 1e-3;
    Def.verbose = 0; Def.plot = 0;
    Def.updateModel = 'new';
%     Def.irka.stopCrit = 's0';
    Def.irka.suppressverbose = 1;
    Def.irka.maxiter = 15;
%     Def.irka.tol = 1e-2; 
    Def.irka.lu = 'full';

    if ~exist('Opts','var') || isempty(Opts)
        Opts = Def;
    else
        Opts = parseOpts(Opts,Def);
    end 
    
    %% run computations
    stop = 0;
    kIter = 0;
    kIrkaTot = 0;
    sysmOld = sss([],[],[]);
    
    %   Generate the model function
    s0m = Opts.s0m;    [sysm, s0mTot, V, W] = modelFct(sys,s0m);

    if Opts.verbose, fprintf('Starting model function MOR...\n'); end

    while ~stop
        kIter = kIter + 1; if Opts.verbose, fprintf(sprintf('modelFctMor: k=%i\n',kIter));end
        if kIter > 1
            % update model
            [sysm, s0mTot, V, W] = modelFct(sys,s0,s0mTot,V,W,Opts);
        end
        % reduction of new model with new starting shifts
        [sysr, ~,~, s0new, ~,~,~,~,~,~,~,kIrka] = irka(sysm,s0,Opts.irka);
        nSysm = sysm.n;
        kIrkaTot = kIrkaTot + kIrka;
        % computation of convergence
        if stoppingCrit
            stop = 1;
        else
            %Overwrite parameters with new variables
            s0 = s0new;    
            sysmOld = sysm;
        end
        if kIter >= Opts.maxiter; 
            warning('modelFctMor did not converge within maxiter'); 
            return
        end
    end

    function stop = stoppingCrit
        stop = 0;
        %   Compute the change in shifts
        normS0 = norm(s0);
        if normS0 == 0
            crit = norm(setdiffVec(s0new,s0)); %absolute
        else
            crit = norm(setdiffVec(s0new,s0))/normS0; %relative
        end
        %   Cumpute the change in model function
        normSysmOld = norm(sysmOld);
        if ~isinf(normSysmOld) %stable
            crit = [crit, norm(sysm-sysmOld)/norm(sysmOld)];
        else
            crit = [crit, NaN];
        end
        
        if any(crit <= Opts.tol), stop = 1;
        elseif length(s0mTot)> size(sys.a,1),stop = 1;end
    end    
end