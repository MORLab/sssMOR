function [sysr, s0, sysm, relH2err] = modelFctMor(sys,redFct,varargin)
% MODELFCTMOR - model function-based model order reduction
%
% Syntax:
%       MODELFCTMOR(sys,redFct)
%       MODELFCTMOR(sys,redFct,s0)
%       MODELFCTMOR(sys,redFct,Opts)
%       MODELFCTMOR(sys,redFct,s0,Opts)
%       sysr = MODELFCTMOR(sys,redFct)
%       [sysr, s0, sysm, relH2err] = MODELFCTMOR(sys,redFct)
%
% Description:
%   This function executes the reduction scheme defined by the handle 
%   redFct, by generating a model function that substitues the original
%   model sys. This substitution saves execution time. 
%   Upon convergence, the model function is updated and the reduction is
%   conducted again until convergence of the overall scheme.
%
%   redFct is a function handle that is given
%       - a dynamical system sys
%       - a set of shifts s0
%   and returns:
%       - a reduced order dynamical system sysr
%       - a new set of (optimal) shifts
%
%   This function is a generalization of |cirka|, which implements the
%   model function framework for H2-optimal reduction using IRKA. In
%   contrast, MODELFCTMOR allows **any** H2-optimal function to be embedded
%   in the model function framework.
%
%   //Note: In its current form, MODELFCTMOR supports only SISO models.
%       An extension to MIMO will be given in a later release.
%
%   //Note: Due to precompilation of functions in newer MATLAB versions,
%   one may experience a worse performance using MODELFCTMOR than e.g.
%   cirka. This is due to the fact that precompilation does not work with
%   function handles.
%
% Input Arguments:  
%       *Required Input Arguments:*
%       -sys:			full oder model (sss)
%       -redFct:        handle to H2-optimal reduction function
%
%       *Optional Input Arguments:*
%       -s0:			vector of initial shifts
%       -Opts:			structure with execution parameters
%			 -.qm0:     initial size of model function;
%                       [{2*length(s0)} / positive integer]
%            -.sm0:     initial shifts for surrogate;
%                       [{[s0,s0]} / vector ]
%            -.maxiter: maximum number of iterations;
%						[{15} / positive integer]
%            -.tol:		convergence tolerance;
%						[{1e-3} / positive float]
%            -.stopcrit:convergence criterion for CIRKA;
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
%
% Output Arguments:      
%       -sysr:              reduced order model (sss)
%       -s0:                final choice of shifts
%       -sysm:              resulting model function
%       -relH2err:          estimate of the relative H2 error
%
% Examples:
%       This code computes an H2-optimal approximation of order 10 to
%       the benchmark model 'building' using Confined IRKA. 
%
%> sys = loadSss('building'); s0 = zeros(1,10);
%> Opts = struct('tol',1e-6);
%> redFct   = @(sys,s) irka(sys,s,Opts);
%> redFctOut= @(sys,s) getDesiredOutput(redFct,[1,4],sys,s);
%> [sysr, s0opt, sysm, relH2err] = modelFctMor(sys,redFctOut,s0);
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
% See also:
%   cirka, irka, cure, spark, modelFct, getDesiredOutput
%
% References:
%       * *[1] Castagnotto et al. (2016)*, Fast H2 optimal model order
%              reduction exploiting the local nature of Krylov-Subspace...
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
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  22 Nov 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

if ~sys.isSiso, error('sssMOR:cirka:notSiso','This function currently works only for SISO models');end
warning('off','Control:analysis:NormInfinite3')

    %%  Input parsing
    %   Varargin
    if nargin > 2
        if ~isstruct(varargin{1})
            s0 = varargin{1};
            if nargin == 4
                Opts = varargin{2};
            end
        else
            s0 = [];
            Opts = varargin{1};
        end
    end

    Def.qm0     = 2*length(s0); %default initial surrogate size
    Def.s0m     = shiftVec([s0;2*ones(1,length(s0))]); %default surrogate shifts
    Def.maxiter = 15;   %maximum number of CIRKA iterations
    Def.tol     = 1e-3; %tolerance for stopping criterion
    Def.stopcrit= 'combAny'; %stopping criterion for CIRKA 
                             %s0,sysr,sysm,combAny,combAll
    Def.verbose = 0; Def.plot = 0; %display text and plots
    Def.suppressWarn = 0; %suppress warnings
    Def.updateModel = 'new'; %shifts used for the model function update
    Def.modelTol = 1e-2; %shift tolerance for model function
    Def.clearInit = 0; %reset the model fct after initialization?
    
    Def.irka.suppressverbose = 1;
    Def.irka.stopcrit        = 'combAny';
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
    if any(strcmp(Opts.stopcrit,{'combAny','combAll'})),
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

    % run model function iteration
    stop = false;
    while ~stop
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
        [sysr, s0new] = redFct(sysm,s0);
        
        if Opts.plot, 
            figure; bodemag(sysFrd,ss(sysm),sysr)
            legend('FOM','ModelFct','ROM');   
            title(sprintf('kIter=%i, nModel=%i',kIter,sysm.n));
            pause
        end
        
        % computation of convergence
        [stop,stopVal(kIter,:)] = stopCrit;
        if Opts.verbose, 
            fprintf(1,'\tstopVal (%s): %s\n',Opts.stopcrit,sprintf('%3.2e\t',stopVal(kIter,:)));
            fprintf(1,'\tModelFct size: %i \n',length(s0mTot));
        end
        
        if stop
            sysm = stabsep(ss(sysm));
            relH2err = norm(ss(sysm-sysr))/norm(ss(sysm));
        else
            %Detect STAGNATION
            if kIter > 1 && max(abs(stopVal(kIter-1,:)-stopVal(kIter,:))) < Opts.tol
                warning('sssMOR:modelFctMor:stagnation','Model Function MOR stagnation: no convergence but no improvement. Aborting.'); 
                sysm = stabsep(ss(sysm));
                relH2err = norm(ss(sysm-sysr))/norm(ss(sysm));
               return
            end
            %Overwrite parameters with new variables
            s0 = s0new;    
            sysmOld = sysm;
            sysrOld = sysr;
        end
        if kIter >= Opts.maxiter; 
            warning('sssMOR:modelFctMor:maxiter','modelFctMor did not converge within maxiter'); 
            if Opts.suppressWarn, warning('on','sssMOR:irka:maxiter');end
            sysm     = stabsep(ss(sysm));
            relH2err = norm(ss(sysm-sysr))/norm(ss(sysm));
            return
        end
    end
    if Opts.suppressWarn,warning('on','sssMOR:irka:maxiter'); end
    if ~Opts.verbose
        fprintf('modelFctMor step %03u - Convergence (%s): %s \n', ...
                kIter, Opts.stopcrit, sprintf('% 3.1e', stopVal(kIter,:)));
    end
    
    % ============== AUXILIARY ===============================
    function [stop,stopVal] = stopCrit(varargin)
        % compute the stopping criterion for CIRKA
        
        % input parsing
        if nargin == 0
            stopcrit = Opts.stopcrit;
        else
            stopcrit = varargin{1};
        end
        
        if length(s0mTot)> sys.n,
            % Full order achieved?
            stop = true; 
        else
            %   Compute the change in shifts
            switch stopcrit
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
                    [stop(1),stopVal(1)] = stopCrit('s0');
                    [stop(2),stopVal(2)] = stopCrit('sysr'); 
                    [stop(3),stopVal(3)] = stopCrit('sysm');
                    stop = all(stop);
                case 'combAny'
                    [stop(1),stopVal(1)] = stopCrit('s0');
                    [stop(2),stopVal(2)] = stopCrit('sysr'); 
                    [stop(3),stopVal(3)] = stopCrit('sysm');
                    stop = any(stop);
                otherwise
                    error('The stopping criterion selected is incorrect or not implemented')
            end
        end
    end    
end