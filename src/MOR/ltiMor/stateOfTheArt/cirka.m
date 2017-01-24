function [sysr, V, W, s0, kIrka, sysm, relH2err] = cirka(sys, s0, Opts) 
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
%       //Note: In its current form, CIRKA supports only SISO models.
%       An extension to MIMO will be given in a later release.
%
% Input Arguments:  
%       *Required Input Arguments:*
%       -sys:			full oder model (sss)
%       -s0:			vector of initial shifts
%
%       *Optional Input Arguments:*
%       -Opts:			structure with execution parameters
%			-.qm0:     initial size of model function;
%                       [{2*length(s0)} / positive integer]
%           -.sm0:     initial shifts for surrogate;
%                       [{[s0,s0]} / vector ]
%           -.maxiter: maximum number of iterations;
%						[{15} / positive integer]
%           -.tol:		convergence tolerance;
%						[{1e-3} / positive float]
%           -.stopcrit:convergence criterion for CIRKA;
%                       ['s0' / 'sysr' / 'sysm' / {'combAny'} / 'combAll']
%           -.verbose:	show text output during iterations;
%						[{false} / true]
%           -.plot:     plot results;
%                       [{false} / true]
%           -.suppressWarn: suppress warnings;
%                       [{false} / true]
%           -.updateModel: type of model function update;
%                       [{'new'},'all']
%           -.clearInit: reset the model function after first iteration;
%                       [{true}, false]
%           -.irka.stopcrit: stopping criterion used in irka;
%                       [{'combAny'} / 's0' / 'sysr' /'combAll']
%           -.irka.lse:  choose type of lse solver;
%                       ['sparse' / {'full'} / 'hess']
%           -.irka.suppressverbose: suppress any type of verbose for speedup;
%                       [{true} / false]
%           -.irka.tol: convergence tolerance for irka
%                       [{1e-6} / positive float ] 
%           (for further irka options, please refer to help irka)
%
% Output Arguments:      
%       -sysr:              reduced order model (sss)
%       -V,W:               resulting projection matrices (V = Vm*Virka)
%       -s0:                final choice of shifts
%       -kIrka:             vector of irka iterations
%       -sysm:              resulting model function
%       -relH2err:          estimate of the relative H2 error
%
% Examples:
%       This code computes an H2-optimal approximation of order 10 to
%       the benchmark model 'building' using Confined IRKA. 
%
%> sys = loadSss('building'); s0 = zeros(1,10);
%> [sysr, ~, ~, s0opt, kIrka, sysm, relH2err] = cirka(sys, s0);
%> bode(sys,'-',sysm,'--r',sysr,'--g'); legend('sys','sysm','sysr')
%
%       //Note: The computational advantage of the model function framework
%       is given especially for truly large scale systems, where the
%       solution of a sparse LSE becomes much more expensive than of a
%       small dense LSE.
%
%       One can use the function |isH2opt| to verify if the necessary 
%       conditions for optimality are satisfied.
%
%> isH2opt(sys, sysr, s0opt);
%
%       Note also that the model function |sysm| allows the estimation of the
%       relative H2 approximation error, given as the last output
%       |relH2err|
%
%> estimatedErr = relH2err
%> actualErr    = norm(sys-sysr)/norm(sys)
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
% Last Change:  20 Jan 2017
% Copyright (c) 2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------
    
%%  Check inputs
if ~sys.isSiso, error('sssMOR:cirka:notSiso','This function currently works only for SISO models');end
warning('off','Control:analysis:NormInfinite3')

%% Define execution options
    Def.qm0     = 2*length(s0); %default initial surrogate size
    Def.s0m     = shiftVec([s0;2*ones(1,length(s0))]); %default surrogate shifts
    Def.maxiter = 15;   %maximum number of CIRKA iterations
    Def.tol     = 1e-3; %tolerance for stopping criterion
    Def.stopCrit= 'combAny'; %stopping criterion for CIRKA 
                             %s0,sysr,sysm,combAny,combAll
    Def.verbose = 0; Def.plot = 0; %display text and plots
    Def.suppressWarn = 0; %suppress warnings
    Def.updateModel = 'new'; %shifts used for the model function update
    Def.modelTol = 1e-2; %shift tolerance for model function
    Def.clearInit = 0; %reset the model fct after initialization?
    
    Def.irka.suppressverbose = true;
    Def.irka.stopCrit        = 'combAny';
    Def.irka.lse             = 'full';
    Def.irka.tol             = 1e-6;

    if ~exist('Opts','var') || isempty(Opts)
        Opts = Def;
    else
        Opts = parseOpts(Opts,Def);
    end 
    
    if Opts.suppressWarn, warning('off','sssMOR:irka:maxiter'); end  
%% run computations
    kIter   = 0;
    kIrka   = zeros(1,Opts.maxiter);
    if any(strcmp(Opts.stopCrit,{'combAny','combAll'})),
        nStopVal = 3;
    else
        nStopVal = 1;
    end
    stopVal = zeros(Opts.maxiter,nStopVal);
    sysrOld = ss([]);
    sysmOld = ss([]);
    
    %   Generate the model function
    s0m = Opts.s0m;    [sysm, s0mTot, Vm, Wm] = modelFct(sys,s0m);

    if Opts.verbose, fprintf('Starting model function MOR...\n'); end
    if Opts.plot, sysFrd = freqresp(sys,struct('frd',true)); end

    % run CIRKA iteration
    stop = false;
    while ~stop &&  kIter <= Opts.maxiter
        kIter = kIter + 1; if Opts.verbose, fprintf(sprintf('modelFctMor: k=%i\n',kIter));end
            
        if kIter > 1
            if kIter == 2 && Opts.clearInit
                %reset the model function after the first step
                s0m = [s0,s0m(1:length(s0m)-length(s0))];
                [sysm, s0mTot, Vm, Wm] = modelFct(sys,s0m);
            else
                % update model
                [sysm, s0mTot, Vm, Wm] = modelFct(sys,s0,s0mTot,Vm,Wm,Opts);
            end
        end
        
        % reduction of new model with new starting shifts
        [sysr, Virka, Wirka, s0new, ~,~,~,~,~,~,~,~,kIrkaNew] = irka(sysm,s0,Opts.irka);

        if Opts.plot, 
            figure; bodemag(sysFrd,ss(sysm),sysr)
            legend('FOM','ModelFct','ROM');   
            title(sprintf('kIter=%i, nModel=%i',kIter,sysm.n));
            pause
        end

        kIrka(kIter) = kIrkaNew;
        
        % computation of convergence
        [stop,stopVal(kIter,:)] = verifyStopCrit;
        %Detect STAGNATION
        if ~stop && kIter > 1 && max(abs(stopVal(kIter-1,:)-stopVal(kIter,:))) < Opts.tol 
            stop = true;
        end
        
        s0      = s0new;    
        sysmOld = sysm;
        sysrOld = sysr;
            
        if Opts.verbose, 
            fprintf(1,'\tkIrka: %03i\n',kIrkaNew);
            fprintf(1,'\tstopVal (%s): %s\n',Opts.stopCrit,sprintf('%3.2e\t',stopVal(kIter,:)));
            fprintf(1,'\tModelFct size: %i \n',length(s0mTot));
        end     
    end
    
    %%  Terminate execution  
    % prepare outputs
    relH2err = norm(sysm-sysr)/norm(sysm);
    kIrka(kIter+1:end) = []; %remove preallocation
    V = Vm*Virka;
    W = Wm*Wirka;
    
    % create ssRed object with CIRKA information
    Opts.s0             = s0;
    Opts.kIrka          = kIrka;
    Opts.originalOrder  = sys.n;
    Opts.modelFctOrder  = sysm.n; 
    
    [Ar,Br,Cr,Dr,Er] = dssdata(sysr);
    sysr = ssRed('cirka',Opts,Ar,Br,Cr,Dr,Er);
    % make model function stable by removing ustable modes
    if ~isstable(sysm)
        n0 = sysm.n;
        sysm = stabsep(sysm);
        % Update ssRed object with information about reduction
        sysm = ssRed('stabsep',struct('originalOrder',n0),sysm);
    end
    
    % display warnings or text output
    if kIter >= Opts.maxiter;
        warning('sssMOR:cirka:maxiter','modelFctMor did not converge within maxiter'); 
    elseif ~Opts.verbose
        fprintf('CIRKA step %03u - Convergence (%s): %s \n', ...
                kIter, Opts.stopCrit, sprintf('% 3.1e', stopVal(kIter,:)));
    end
    
    if Opts.suppressWarn,warning('on','sssMOR:irka:maxiter'); end

     
% ============== AUXILIARY ===============================
    function [stop,stopVal] = verifyStopCrit(varargin)
        % compute the stopping criterion for CIRKA
        
        % input parsing
        if nargin == 0
            stopCrit = Opts.stopCrit;
        else
            stopCrit = varargin{1};
        end
        
        if length(s0mTot)> sys.n,
            % Full order achieved?
            stop = true; 
        else
            %   Compute the change in shifts
            switch stopCrit
                case 's0' 
                    %shift convergence
                    if any(abs(s0))<1e-3
                        stopVal = norm((s0-s0new), 1)/sysr.n;
                    else
                        stopVal = norm((s0-s0new)./s0, 1)/sysr.n;
                    end
                    stop = (stopVal <= Opts.tol);
                    
                case 'sysr' 
                    %reduced model convergence
%                     stopVal = inf; %initialize in case the reduced model is unstable
%                     if all(real(eig(sysr))<0) && all(real(eig(sysrOld))<0)
                    stopVal=norm(sysr-sysrOld)/norm(sysr);
%                     end
                    stop = (stopVal <= Opts.tol);
                case 'sysm' 
                    %model function convergence
                    stopVal=norm(sysm-sysmOld)/norm(sysm);
                    stop = (stopVal <= Opts.tol);
                case 'combAll'
                    [stop(1),stopVal(1)] = verifyStopCrit('s0');
                    [stop(2),stopVal(2)] = verifyStopCrit('sysr'); 
                    [stop(3),stopVal(3)] = verifyStopCrit('sysm');
                    stop = all(stop);
                case 'combAny'
                    [stop(1),stopVal(1)] = verifyStopCrit('s0');
                    [stop(2),stopVal(2)] = verifyStopCrit('sysr'); 
                    [stop(3),stopVal(3)] = verifyStopCrit('sysm');
                    stop = any(stop);
                otherwise
                    error('The stopping criterion selected is incorrect or not implemented')
            end
        end
    end    
end
