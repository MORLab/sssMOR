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
            Opts = varargin{1};
        end
    end

    %   Default Opts
    Def.qm0     = 10;
    Def.s0m     = zeros(1,Def.qm0);
    Def.maxiter = 20;
    Def.tol     = 1e-6;
    Def.verbose = 0;

    if ~exist('Opts','var') || isempty(Opts)
        Opts = Def;
    else
        Opts = parseOpts(Opts,Def);
    end  
    %%  Computations

    s0m = Opts.s0m; %initialize
    stop = 0;
    kIter = 0;

    if Opts.verbose, fprintf('Starting model function MOR...\n'); end

    while ~stop
        kIter = kIter + 1; if Opts.verbose, fprintf(sprintf('k=%i',k));end
        % update model
        sysm = rk(sys,s0m,s0m);
        % reduction of new model with new starting shifts
        [sysr, s0new] = redFct(sysm,s0);
        % computation of convergence
        if stoppingCrit
            stop = 1;
        else
            s0m = [s0m,s0new]; %shifts for updated model
            s0 = s0new;
        end
        if kIter > Opts.maxiter; 
            error('modelFctMor did not converge within maxiter'); 
        end
    end



    function stop = stoppingCrit
        stop = 0;
        if norm(setdiffVec(s0new,s0))/norm(s0)<= Opts.tol, stop = 1; end
    end
end