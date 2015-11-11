function [sysr, V, W, s0, s0Traj, Rsylv, Lsylv] = irka(sys, s0, varargin) 
% IRKA - Iterative Rational Krylov Algorithm
%
% Syntax:
%       sysr = IRKA(sys, s0)
%       sysr = IRKA(sys, s0, Opts)
%       sysr = IRKA(sys, s0, Rt, Lt)
%       sysr = IRKA(sys, s0, Rt, Lt, Opts)
%       [sysr, V, W] = IRKA(sys, ... )
%       [sysr, V, W, s0] = IRKA(sys, ... )
%       [sysr, V, W, s0, s0Traj] = IRKA(sys, ... )
%       [sysr, V, W, s0, s0Traj, Rsylv, Lsylv] = IRKA(sys, ... )
%
% Description:
%       This function executes the Iterative Rational Krylov
%       Algorithm (IRKA) as proposed by Gugergin and Beattie in [1].
% 
%       The IRKA iteration is conducted to search for an optimal set of 
%       shifts in Krylov subspace-based model reduction. If IRKA converges,
%       then the reduced model is known to be a local optimum with respect
%       to the H2 norm of the error.
%
% Input Arguments:       
%       -sys:       full oder model (sss)
%       -s0:        vector of initial shifts
%       -Opts:      (opt.) structure with execution parameters
%       -Rt/Lt:     initial right/left tangential directions for MIMO
%
% Output Arguments:      
%       -sysr:     reduced order model (sss)
%       -V,W:      resulting projection matrices
%       -s0:       final choice of shifts
%       -s0Traj:  trajectory of all shifst for all iterations
%       -Rt,Lt:   right/left tangential directions from arnoldi
%
% Examples:
%       TODO
%
% See Also: 
%       arnoldi, rk
%
% References:
%       * *[1] Gugercin (2008)*, H2 model reduction for large-scale linear dynamical systems
%       * *[2] Beattie (2014)*, Model reduction by rational interpolation
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
% Authors:      Heiko Panzer, Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  31 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Parse input and load default parameters
if nargin > 2
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
    elseif nargin == 5
        %usage irka(sys,s0,Rt,Lt,Opts)
        Rt = varargin{1};
        Lt = varargin{2};
        Opts = varargin{3};
    end
end
if sys.isMimo && (~exist('Rt','var')|| isempty(Rt) ...
        || ~exist('Lt','var')|| isempty(Lt))
    error('specify initial tangential directions for MIMO systems');
end
    
%% Parse the inputs
%   Default execution parameters
Def.maxiter = 50; 
Def.tol = 1e-3; 
Def.type = ''; %'stab', 'newton', 'restarted'
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
s0Traj = zeros(Opts.maxiter+2, length(s0));
s0Traj(1,:) = s0;

%% IRKA iteration
k=0;
while true
    k=k+1; sysr_old = sysr;
    %   Reduction
    if sys.isSiso
        [sysr, V, W, ~ , Rsylv, ~ , Lsylv] = rk(sys, s0, s0);
    else
        [sysr, V, W, ~, Rsylv, ~ , Lsylv] = rk(sys, s0, s0, Rt, Lt);
    end 
    %--
%     warning off
%     rMatching = zeros(1,sysr.n); lMatching = zeros(1,sysr.n);
%     for iS = 1:length(s0)
%         rMatching(iS) = (norm(moments(sys-sysr,s0(iS),1)*Rt(:,iS)));
%         lMatching(iS) = (norm(Lt(:,iS).'*moments(sys-sysr,s0(iS),1)));
%     end
%     all([rMatching, lMatching]< 1e-5)
%     for iS = 1:length(s0)
%         rMatching(iS) = (norm(moments(sys-sysr,s0(iS),1)*Rsylv(:,iS)));
%         lMatching(iS) = (norm(Lsylv(:,iS).'*moments(sys-sysr,s0(iS),1)));
%     end
%     all([rMatching, lMatching]< 1e-5)
%     warning on
    %--
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
    s0Traj(k+1,:) = s0;
    
    [stop, stopCrit] = stoppingCriterion(s0,s0_old,sysr,sysr_old,Opts);
    if Opts.verbose
        fprintf('IRKA step %03u - Convergence: %s \n', ...
            k, sprintf('% 3.1e', stopCrit));
    end
    if stop || k>= Opts.maxiter
        s0 = s0_old; % function return value
        s0Traj = s0Traj(1:(k+1),:);
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

    
    




