function [sysr, V, W, s0, s0Traj, Rt, Lt, B_, SRsylv, Rsylv, C_, SLsylv, Lsylv, kIter] = irka(sys, s0, varargin) 
% IRKA - Iterative Rational Krylov Algorithm
%
% Syntax:
%       sysr                            = IRKA(sys, s0)
%       sysr                            = IRKA(sys, s0)
%       sysr                            = IRKA(sys, s0, Rt, Lt)
%       sysr                            = IRKA(sys, s0,..., Opts)
%       [sysr, V, W]                    = IRKA(sys, s0,... )
%       [sysr, V, W, s0]                = IRKA(sys, s0,... )
%       [sysr, V, W, s0, s0Traj]        = IRKA(sys, s0,... )
%       [sysr, V, W, s0, s0Traj, Rt, Lt]= IRKA(sys, s0,... )
%       [sysr, V, W, s0, s0Traj, Rt, Lt, B_, SRslyv, Rsylv, C_, SLsylv, Lsylv, kIter] = IRKA(sys, s0,... )
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
%       Convergence is determined by observing the shifts and norm of the
%       reduced model over the iterations. This behavior can be changed
%       with the optional Opts structure.
%
% Input Arguments:  
%       *Required Input Arguments:*
%       -sys:			full oder model (sss)
%       -s0:			vector of initial shifts
%
%       *Optional Input Arguments:*
%       -Rt/Lt:			initial right/left tangential directions for MIMO
%       -Opts:			structure with execution parameters
%			-.maxiter:	maximum number of iterations;
%						[{50} / positive integer]
%			-.tol:		convergence tolerance;
%						[{1e-3} / positive float]
%			-.type:		choose between different irka modifications;
%						[{''} / 'stab']
%			-.verbose:	show text output during iterations;
%						[{0} / 1]
%			-.stopCrit:	stopping criterion;
%						[{'combAny'} / 's0' / 'sysr' / 'combAll']
%           -.suppressverbose: suppress any type of verbose for speedup;
%                       [{0} / 1]
%           -.real:     keep the projection matrices real
%                       [{'real'} / '0']
%           -.orth:     orthogonalization of new projection direction
%                       [{'2mgs'} / 0 / 'dgks' / 'mgs']
%           -.reorth:   reorthogonalization
%                       [{'gs'} / 0 / 'qr']
%           -.lse:      use LU or hessenberg decomposition
%                       [{'sparse'} / 'full' / 'hess']
%           -.dgksTol:  tolerance for dgks orthogonalization
%                       [{1e-12} / positive float]
%           -.krylov:   standard or cascaded krylov basis
%                       [{0} / 'cascade]
%
% Output Arguments:      
%       -sysr:              reduced order model (sss)
%       -V,W:               resulting projection matrices
%       -s0:                final choice of shifts
%       -s0Traj:            trajectory of all shifst for all iterations
%       -B_,SRsylv,Rsylv:   matrices of the input Sylvester equation
%       -C_,SLsylv,Lsylv:   matrices of the output Sylvester equation
%       -kIter:             number of iterations
%
% Examples:
%       This code computes an H2-optimal approximation of order 8 to
%       the benchmark model 'fom'. One can use the function isH2opt to
%       verify if the necessary conditions for optimality are satisfied.
%
%> sys = loadSss('fom')
%> [sysr, ~, ~, s0opt] = irka(sys, -eigs(sys,8).');
%> bode(sys,'-',sysr,'--r');
%> isH2opt(sys, sysr, s0opt)
%
% See Also: 
%       arnoldi, rk, isH2opt
%
% References:
%       * *[1] Gugercin et al. (2008)*, H2 model reduction for large-scale linear dynamical systems
%       * *[2] Beattie et al. (2014)*, Model reduction by rational interpolation
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
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin = varargin(1:end-1);
else
    Opts = struct();
end
if ~isempty(varargin)
        %usage irka(sys,s0,Rt,Lt)
        Rt = varargin{1};
        Lt = varargin{2};
elseif ~sys.isSiso
        error('specify initial tangential directions for MIMO systems');
end
    
%% Parse the inputs
%   Default execution parameters
Def.maxiter = 50; 
Def.tol = 1e-3; 
Def.type = ''; %'stab', 'newton', 'restarted'
Def.verbose = 0; % text output durint iteration?
Def.stopCrit = 'combAny'; %'s0', 'sysr', 'combAll', 'combAny'
Def.suppressverbose = 0;

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

% sort expansion points & tangential directions
s0old = s0;
s0 = cplxpair(s0);
if exist('Rt','var') && ~isempty(Rt)
    [~,cplxSorting] = ismember(s0,s0old); 
    Rt = Rt(:,cplxSorting);
    Lt = Lt(:,cplxSorting);
else
    Rt = ones(sys.m,length(s0));
    Lt = ones(sys.p,length(s0));
end
clear s0old

% Initialize variables
sysr = sss([],[],[]);
s0Traj = zeros(Opts.maxiter+2, length(s0));
s0Traj(1,:) = s0;

%% IRKA iteration
kIter=0;
while true
    kIter=kIter+1; sysr_old = sysr;
    %   Reduction
    if sys.isSiso
        [sysr, V, W, B_, SRsylv, Rsylv, C_, SLsylv, Lsylv] = rk(sys, s0, s0,Opts);
    else
        [sysr, V, W, B_, SRsylv, Rsylv, C_, SLsylv, Lsylv] = rk(sys, s0, s0, Rt, Lt,Opts);
    end 
    
    %   Update of the reduction parameters
    s0_old=s0; if ~sys.isSiso, Rt_old = Rt; Lt_old = Lt; end
    if sys.isMimo
        [X, D, Y] = eig(sysr);
        Rt = full((Y.'*sysr.B).'); Lt = full(sysr.C*X);
        s0 = full(-diag(D).');
        % make sure real shifts have real directions
        idx = find(imag(s0)==0);
        if any([max(imag(Rt(:,idx))),max(imag(Lt(:,idx)))]) > 1e-15
            warning('Tangential directions corresponding to real shifts are complex')              
        end
        Rt(:,idx) = real(Rt(:,idx)); Lt(:,idx) = real(Lt(:,idx));
    else
        s0 = -eig(sysr)';
    end
    if strcmp(Opts.type,'stab')
        % mirror shifts with negative real part
        s0 = s0.*sign(real(s0));
    end
    s0Traj(kIter+1,:) = s0;
    
    [stop, stopCrit] = stoppingCriterion(s0,s0_old,sysr,sysr_old,Opts);
    if Opts.verbose
        fprintf('IRKA step %03u - Convergence: %s \n', ...
            kIter, sprintf('% 3.1e', stopCrit));
    end
    if stop || kIter>= Opts.maxiter
        s0 = s0_old; % function return value
        if ~sys.isSiso, Rt = Rt_old; Lt = Lt_old; end
        s0Traj = s0Traj(1:(kIter+1),:);
        sysr.Name = sprintf('%s_irka',sys.Name);
        break
    end      
end
if ~Opts.suppressverbose %display at least the last value
    fprintf('IRKA step %03u - Convergence (%s): %s \n', ...
            kIter, Opts.stopCrit, sprintf('% 3.1e', stopCrit));
end
if kIter>=Opts.maxiter
    warning('sssMOR:irka:maxiter',['IRKA has not converged after ' num2str(kIter) ' steps.']);
    return
end

%%------------------ AUXILIARY FUNCTIONS -------------------
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

    
    




