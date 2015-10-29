function [sysr, V, W, s0, s0_traj] = irka(sys, s0, varargin) 
% IRKA - Iterative Rational Krylov Algorithm
% ------------------------------------------------------------------
% [sysr, V, W, s0, s0_traj] = IRKA(sys, s0, Opts)
% Inputs:       * sys: full oder model (sss)
%               * s0:  vector of initial shifts
%               * Opts (opt.) structure with execution parameters
% Outputs:      * sysr: reduced order model (sss)
%               * V,W:  resulting projection matrices
%               * s0:   final choice of shifts
%               * s0_traj: trajectory of all shifst for all iterations
% ------------------------------------------------------------------
% USAGE:  This function executes the Iterative Rational Krylov
% Algorithm (IRKA) as proposed by Gugergin and Beattie in [1].
%
% The IRKA iteration is conducted to search for an optimal set of 
% shifts in Krylov subspace-based model reduction. If IRKA converges,
% then the reduced model is known to be a local optimum with respect
% to the H2 norm of the error.
%
% See also ARNOLDI, RK.
%
% ------------------------------------------------------------------
% REFERENCES:
% [1] Gugercin (2008), H2 model reduction for large-scale linear
%     dynamical systems
% [2] Beattie (2014), Model reduction by rational interpolation
% ------------------------------------------------------------------
% This file is part of MORLab, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute 
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                   -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Heiko Panzer, Alessandro Castagnotto
% Last Change:  28 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%% Parse input and load default parameters
if nargin == 3
    %usage irka(sys,s0,Opts)
    Opts = varargin{1};
    if sys.isMimo
        error('specify initial tangential directions for MIMO systems');
    end
elseif nargin == 4
    %usage irka(sys,s0,Rt,Lt)
    Rt = varargin{1};
    Lt = varargin{2};
else
    %usage irka(sys,s0,Rt,Lt,Opts)
    Rt = varargin{1};
    Lt = varargin{2};
    Opts = varargin{1};
end
%% Parse the inputs
%   Default execution parameters
Def.maxiter = 50; 
Def.tol = 1e-3; 
Def.type = ''; %'stab', 'newton','restarted'
Def.verbose = 0; % text output durint iteration?
Def.stopCrit = 'combAny'; %'s0', 'sysr', 'combAll', 'combAny'
Def.cplxpairTol = 1e-6;

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

% Further check/transform the input
if Opts.tol<=0 || ~isreal(Opts.tol)
    error('tol must be a real positive number.');
end

s0 = s0_vect(s0);

% Initialize variables
sysr = sss([],[],[]);
s0_traj = zeros(Opts.maxiter+2, length(s0));
s0_traj(1,:) = s0;

%% IRKA iteration
k=0;
while true
    k=k+1; sysr_old = sysr;
    %   Reduction
    if sys.isMimo
        [sysr, V, W] = rk(sys, s0, s0, Rt, Lt);
    else
        [sysr, V, W] = rk(sys, s0, s0);
    end
    %   Update of the reduction parameters
    s0_old=s0;
    if sys.isMimo
        [X, D, Y] = eig(sysr);
        Rt = full((Y.'*sysr.B).'); Lt = full(sysr.C*X);
        s0 = full(-diag(D).');
    else
        s0 = -eig(sysr)';
    end
    if strcmp(Opts.type,'stab')
        % mirror shifts with negative real part
        s0 = s0.*sign(real(s0));
    end
    s0_traj(k+1,:) = s0;
    
    [stop, stopCrit] = stoppingCriterion(s0,s0_old,sysr,sysr_old,Opts);
    if Opts.verbose
        fprintf('IRKA step %03u - Convergence: %s \n', ...
            k, sprintf('% 3.1e', stopCrit));
    end
    if stop || k>= Opts.maxiter
        s0 = s0_old; % function return value
        s0_traj = s0_traj(1:(k+1),:);
        break
    end      
end
if ~Opts.verbose %display at least the last value
    fprintf('IRKA step %03u - Convergence (%s): %s \n', ...
            k, Opts.stopCrit, sprintf('% 3.1e', stopCrit));
end
if k>=Opts.maxiter
    warning('IRKA:no_converged', ['IRKA has not converged after ' num2str(k) ' steps.']);
    return
end

%------------------ AUXILIARY FUNCTIONS -------------------
function s0=s0_vect(s0)
    % change two-row notation to vector notation
    if size(s0,1)==2
        temp=zeros(1,sum(s0(2,:)));
        for j=1:size(s0,2)
            k=sum(s0(2,1:(j-1))); k=(k+1):(k+s0(2,j));
            temp(k)=s0(1,j)*ones(1,s0(2,j));
        end
        s0=temp;
    end
    % sort expansion points
    s0 = cplxpair(s0);
    if size(s0,1)>size(s0,2)
        s0=transpose(s0);
    end
function [stop,stopCrit] = stoppingCriterion(s0,s0_old,sysr,sysr_old,Opts)
%   Computes the stopping criterion(s) for IRKA
%   It returns a scalar or a vector depending on the number of criteria
%   chosen
switch Opts.stopCrit
    case 's0' %shift convergence
        stopCrit = norm((s0-s0_old)./s0, 1)/sysr.n;
        stop = stopCrit <= Opts.tol;
    case 'sysr' %reduced model convergence
        stopCrit = inf; %initialize in case the reduced model is unstable
        if all(real(eig(sysr))<0) && all(real(eig(sysr_old))<0)
                stopCrit=norm(sysr-sysr_old)/norm(sysr);
        end
        stop = stopCrit <= Opts.tol;
    case 'combAll'
        stopCrit = norm((s0-s0_old)./s0, 1)/sysr.n;
        stopCrit = [stopCrit, inf]; 
        if all(real(eig(sysr))<0) && all(real(eig(sysr_old))<0)
                stopCrit(2) = norm(sysr-sysr_old)/norm(sysr);
        end
        stop = all(stopCrit <= Opts.tol);
    case 'combAny'
        stopCrit = norm((s0-s0_old)./s0, 1)/sysr.n;
        stopCrit = [stopCrit, inf]; 
        if all(real(eig(sysr))<0) && all(real(eig(sysr_old))<0)
                stopCrit(2) = norm(sysr-sysr_old)/norm(sysr);
        end
        stop = any(stopCrit <= Opts.tol);
    otherwise
        error('The stopping criterion selected is incorrect or not implemented')
end

    
    




