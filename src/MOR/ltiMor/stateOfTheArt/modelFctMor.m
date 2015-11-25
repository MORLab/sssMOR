function [sysr, s0new] = modelFctMor(sys,redFct,varargin)
% MODELFCTMOR - model function-based model order reduction
%
% Syntax:
%   MODELFCTMOR(sys,redFct)
%   MODELFCTMOR(sys,redFct,s0)
%   MODELFCTMOR(sys,redFct,Opts)
%   MODELFCTMOR(sys,redFct,s0,Opts)
%   sysr = MODELFCTMOR(sys,redFct)
%   [sysr, s0] = MODELFCTMOR(sys,redFct)
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
% See also:
%   IRKA, CURE, SPARK
%
% References:
%   [1] Castagnotto (2015), tbd
%
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Last Change:  01 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

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

    %   Default Opts
    Def.qm0     = length(s0)+2; %at least reduced order + 2
    Def.s0m     = zeros(1,Def.qm0);
    Def.maxiter = 8;
    Def.tol     = 1e-3;
    Def.verbose = 0;
    Def.updateModel = 'new'; % 'all','new','lean'

    if ~exist('Opts','var') || isempty(Opts)
        Opts = Def;
    else
        Opts = parseOpts(Opts,Def);
    end  
    %%  Computations
    
    %initialize variables
    stop = 0;
    kIter = 0;
    
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
        [sysr, s0new] = redFct(sysm,s0);
        % computation of convergence
        if stoppingCrit
            stop = 1;
        else
            %Overwrite parameters with new variables
            s0 = s0new;    
        end
        if kIter > Opts.maxiter; 
            warning('modelFctMor did not converge within maxiter'); 
            return
        end
    end

    function stop = stoppingCrit
        stop = 0;
        if norm(s0) == 0
            crit = norm(setdiffVec(s0new,s0)); %absolute
        else
            crit = norm(setdiffVec(s0new,s0))/norm(s0); %relative
        end
        if crit <= Opts.tol, stop = 1;
        elseif length(s0mTot)> size(sys.a,1),stop = 1;end
    end    
end