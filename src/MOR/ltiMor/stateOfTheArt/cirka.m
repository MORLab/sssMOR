function [sysr, Virka, Wirka, s0, kIrka, sysm, relH2err] = cirka(sys, s0, Opts) 
% CIRKA - Confined Iterative Rational Krylov Algorithm
%
% Syntax:
%       sysr                                    = CIRKA(sys, s0)
%       sysr                                    = CIRKA(sys, s0)
%       sysr                                    = CIRKA(sys, s0, Opts)
%       [sysr, V, W]                            = CIRKA(sys, s0,... )
%       [sysr, V, W, s0]                        = CIRKA(sys, s0,... )
%       [sysr, V, W, s0, kIrka, sysm, relH2err] = CIRKA(sys, s0,... )
%
% Description:
%       This function executes the Confined Iterative Rational Krylov
%       Algorithm (CIRKA) as proposed by Castagnotto et al. in [1].
% 
%       The algorithm is based on constructing a model function, i.e. a
%       surrogate representing the full oder model locally about some
%       frequencies, and running IRKA [2] with respect to the surrogate model.
%       The model function is updated until convergence.
%
%       Convergence is determined by observing the shifts and norm of the
%       reduced model over the iterations. This behavior can be changed
%       with the optional Opts structure.
%
%       *In its current form, CIRKA supports only SISO models*
%       An extension will be given in a later release.
%
% Input Arguments:  
%       *Required Input Arguments:*
%       -sys:			full oder model (sss)
%       -s0:			vector of initial shifts
%
%       *Optional Input Arguments:*
%       -Opts:			structure with execution parameters
%			 -.qm0     = initial size of model function;
%                       [{2*length(s0)} / positive integer]
%            -.sm0     = initial shifts for surrogate;
%                       [{[s0,s0]} / vector ]
%            -.maxiter = maximum number of iterations;
%						[{15} / positive integer]
%            -.tol:		convergence tolerance;
%						[{1e-3} / positive float]
%           -.verbose:	show text output during iterations;
%						[{false} / true]
%           -.plot:     plot results;
%                       [{false} / true]
%           -.suppressWarn: suppress warnings;
%                       [{false} / true]
%           -.updateModel: type of model function update;
%                       [{'new'},'all']
%           -.clearInit: reset the model function after first iteration
%                       [{true}, false]
%           -irka.stopcrit: stopping criterion used in irka
%                       [{'combAny'} / 's0' / 'sysr' /'combAll']
%           -irka.lse:  choose type of lse solver
%                       ['sparse' / {'full'} / 'hess']
%           -irka.suppressverbose: suppress any type of verbose for speedup;
%                       [{0} / 1]
%           (for further irka options, please refer to help irka)
%
% Output Arguments:      
%       -sysr:              reduced order model (sss)
%       -V,W:               resulting projection matrices
%       -s0:                final choice of shifts
%       -kIrka:             vector of irka iterations
%       -sysm:              resulting model function
%       -relH2err:          estimate of the relative H2 error
%
% Examples:
%       This code computes an H2-optimal approximation of order 8 to
%       the benchmark model 'fom'. One can use the function isH2opt to
%       verify if the necessary conditions for optimality are satisfied.
%
%> sys = loadSss('fom')
%> [sysr, ~, ~, s0opt] = cirka(sys, -eigs(sys,8).');
%> bode(sys,'-',sysr,'--r');
%> isH2opt(sys, sysr, s0opt)
%
% See Also: 
%       arnoldi, rk, isH2opt, irka, modelFct
%
% References:
%       * *[1] Castagnotto et al. (2016)*, Fast H2 optimal model order
%              reduction exploiting the local nature of Krylov-Subspace...
%       * *[2] Gugercin et al. (2008)*, H2 model reduction for large-scale linear dynamical systems
%       * *[3] Beattie et al. (2014)*, Model reduction by rational interpolation
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  20 Nov 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------
    
if ~sys.isSiso, error('sssMOR:cirka:notSiso','This function currently works only for SISO models');end

%% Define execution options
    Def.qm0     = 2*length(s0); %default initial surrogate size
    Def.s0m     = shiftVec([s0;2*ones(1,length(s0))]); %default surrogate shifts
    Def.maxiter = 15;   %maximum number of CIRKA iterations
    Def.tol     = 1e-3; %tolerance for stopping criterion
    Def.verbose = 0; Def.plot = 0; %display text and plots
    Def.suppressWarn = 0; %suppress warnings
    Def.updateModel = 'new'; %shifts used for the model function update
    Def.modelTol = 1e-1; %shift tolerance for model function
    Def.clearInit = 0; %reset the model fct after initialization?
    
    Def.irka.suppressverbose = 1;
    Def.irka.stopcrit        = 'combAny';
    Def.irka.lse             = 'full';

    if ~exist('Opts','var') || isempty(Opts)
        Opts = Def;
    else
        Opts = parseOpts(Opts,Def);
    end 
    
    if Opts.suppressWarn, warning('off','sssMOR:irka:maxiter'); end  
%% run computations
    stop = 0;
    kIter = 0;
    kIrka = [];
%     sysmOld = sss([],[],[]);
    
    %   Generate the model function
    s0m = Opts.s0m;    [sysm, s0mTot, V, W] = modelFct(sys,s0m);

    if Opts.verbose, fprintf('Starting model function MOR...\n'); end
    if Opts.plot, sysFrd = freqresp(sys,struct('frd',true)); end

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
        [sysr, Virka, Wirka, s0new, ~,~,~,~,~,~,~,~,kIrkaNew] = irka(sysm,s0,Opts.irka);

        if Opts.plot, 
            fh = figure; bodemag(sysFrd,ss(sysm),sysr)
            legend('FOM','ModelFct','ROM');   
            title(sprintf('kIter=%i, nModel=%i',kIter,sysm.n));
            pause
        end
%         if Opts.plot, 
%             figure(fh); plot(real(s0),imag(s0),'bx'); hold on
%                        plot(real(s0new),imag(s0new),'or'); hold off
%                        xlabel('Re'); ylabel('Im'); title('CIRKA shifts')
%                        legend('old','new')
%                        pause
%         end

        kIrka = [kIrka, kIrkaNew];
        % computation of convergence
        if stoppingCrit
            stop = true;
            sysm = stabsep(ss(sysm));
            relH2err = norm(ss(sysm-sysr))/norm(ss(sysm));
        else
            %Overwrite parameters with new variables
            s0 = s0new;    
            sysmOld = sysm;
        end
        if kIter >= Opts.maxiter; 
            warning('modelFctMor did not converge within maxiter'); 
            if Opts.suppressWarn, warning('on','sssMOR:irka:maxiter');end
            sysm = stabsep(ss(sysm));
            relH2err = norm(ss(sysm-sysr))/norm(ss(sysm));
            return
        end
    end
    if Opts.suppressWarn,warning('on','sssMOR:irka:maxiter'); end
    if ~Opts.verbose
    fprintf('CIRKA step %03u - Convergence: %s \n', ...
            kIter, sprintf('% 3.1e', crit(1)));
    end
       
    function stop = stoppingCrit
        stop = false;
        %   Compute the change in shifts
        if any(abs(s0))<1e-3
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
            crit = crit;
%         end
%         crit = [crit, inf]; %use only s0 for convergence
        if Opts.verbose, 
            fprintf(1,'\tcrit: %f\n',crit(1));
            fprintf(1,'\tModelFct size: %i \n',length(s0mTot));
        end
        
        if any(crit <= Opts.tol), stop = 1;
        % Full order achieved?
        elseif length(s0mTot)> sys.n,stop = 1; end
    end    
end
