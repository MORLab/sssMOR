function [sysr, s0, kIter, kIrkaTot, nSysm] = cirka(sys, s0, Opts) 

    %% Define execution options
    Def.qm0     = length(s0)+2;
    Def.s0m     = zeros(1,Def.qm0); 
    Def.maxiter = 15; Def.tol = 1e-3;
    Def.verbose = 0; Def.plot = 0;
    Def.suppressWarn = 0;
    Def.updateModel = 'new';
    Def.clearInit = 1; %reset the model fct after initialization
    Def.irka.stopCrit = 's0';
    Def.irka.suppressverbose = 1;
    Def.irka.maxiter = 100;
%     Def.irka.tol = 1e-2; 
    Def.irka.lse = 'full';

    if ~exist('Opts','var') || isempty(Opts)
        Opts = Def;
    else
        Opts = parseOpts(Opts,Def);
    end 
    
    if Opts.suppressWarn, warning('off','sssMOR:irka:maxiter'); end
    
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
            if kIter == 2 && Opts.clearInit
                %reset the model function after the first step
                s0m = [s0,s0m(1:length(s0m)-length(s0))];
                [sysm, s0mTot, V, W] = modelFct(sys,s0m);
            else
                % update model
                [sysm, s0mTot, V, W] = modelFct(sys,s0,s0mTot,V,W,Opts);
            end
        end
        % reduction of new model with new starting shifts
        [sysr, ~,~, s0new, ~,~,~,~,~,~,~,~,~,kIrka] = irka(sysm,s0,Opts.irka);
        if Opts.plot, 
            figure(1); plot(real(s0),imag(s0),'bx'); hold on
                       plot(real(s0new),imag(s0new),'or'); hold off
                       pause
        end
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
            if Opts.suppressWarn, warning('on','sssMOR:irka:maxiter');end
            return
        end
    end
    if Opts.suppressWarn,warning('on','sssMOR:irka:maxiter'); end
    if ~Opts.verbose
    fprintf('CIRKA step %03u - Convergence: %s \n', ...
            kIter, sprintf('% 3.1e', crit(1)));
    end

        
    function stop = stoppingCrit
        stop = 0;
        %   Compute the change in shifts
        normS0 = norm(s0);
        if normS0 == 0
%             crit = norm(setdiffVec(s0new,s0)); %absolute
            crit = norm((s0-s0new), 1)/sysr.n;
        else
%             crit = norm(setdiffVec(s0new,s0))/normS0; %relative
            crit = norm((s0-s0new)./s0, 1)/sysr.n;
        end
        %   Compute the change in model function
%         normSysmOld = norm(sysmOld);
%         if ~isinf(normSysmOld) %stable
%             crit = [crit, norm(sysm-sysmOld)/norm(sysmOld)];
%         else
            crit = [crit, NaN];
%         end
%         crit = [crit, inf]; %use only s0 for convergence
        if Opts.verbose, 
            fprintf(1,'\tcrit: [%f, %f]\n',crit(1),crit(2));
            fprintf(1,'\tModelFct size: %i \n',length(s0mTot));
        end
        
        if any(crit <= [Opts.tol, 1e-9]), stop = 1;
        % Full order achieved?
        elseif length(s0mTot)> sys.n,stop = 1; end
    end    
end
