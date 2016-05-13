function [sysr, V, W, s0, Rt, B_, Rsylv, s0Traj, RtTraj] = isrk(sys, s0, varargin) 
% ISRK - Iterative SVD-Rational Krylov Algorithm
%
% Syntax:
%       sysr                            = ISRK(sys, s0)
%       sysr                            = ISRK(sys, s0)
%       sysr                            = ISRK(sys, s0, Rt)
%       sysr                            = ISRK(sys, s0,..., Opts)
%       [sysr, V, W]                    = ISRK(sys, s0,... )
%       [sysr, V, W, s0]                = ISRK(sys, s0,... )
%       [sysr, V, W, s0, Rt,]           = ISRK(sys, s0,... )
%       [sysr, V, W, s0, Rt, B_, Rsylv] = ISRK(sys, s0,... )
%       [sysr, V, W, s0, Rt, B_, Rsylv, s0Traj, RtTraj] = ISRK(sys, s0,... )
%
%
% Description:
%       This function executes the Iterative SVD-Rational Krylov
%       Algorithm (ISRK) as proposed by Gugergin in [1].
% 
%       The ISRK iteration combines SVD-based methods like balanced
%       truncation and rational Krylov. The algorithm requires computing
%       the observability Gramian and running an iteration to find a set of
%       pseudo-optimal shifts.
%
%       If ISRK converges, then the reduced model is known to be a local 
%       pseudo-optimum with respect to the H2 norm of the error,
%       pseudo-optimum meaning that it is optimal in the subset of reduced
%       order models with same spectrum and input residues. 
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
%       -Rt:			initial right/left tangential directions for MIMO
%       -Opts:			structure with execution parameters
%			-.maxiter:	maximum number of iterations;
%						[{50} / positive integer]
%			-.tol:		convergence tolerange;
%						[{1e-3} / positive float]
%			-.verbose:	show text output during iterations;
%						[{0} / 1]
%			-.stopCrit:	stopping criterion;
%						[{'combAny'} / 's0' / 'sysr' / 'combAll']
%
% Output Arguments:      
%       -sysr:			reduced order model (sss)
%       -V,W:			resulting projection matrices
%       -s0:			final choice of shifts
%       -s0Traj:		trajectory of all shifst for all iterations
%
% Examples:
%       This code computes an H2-pseudo-optimal approximation of order 8 to
%       the benchmark model 'fom'. 
%
%> sys = loadSss('fom')
%> [sysr, ~, ~, s0opt] = isrk(sys, -eigs(sys,8).');
%> bode(sys,'-',sysr,'--r');
%
% See Also: 
%       arnoldi, rk, irka, tbr
%
% References:
%       * *[1] Gugercin (2008)*, An iterative SVD-Krylov based method for
%              model reduction of large-scale dynamical systems
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
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  23 Dec 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------


%% Parse the inputs

if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin = varargin(1:end-1);
else
    Opts = struct();
end
if ~isempty(varargin)
        %usage isrk(sys,s0,Rt)
        Rt = varargin{1};
elseif ~sys.isSiso
        error('specify initial tangential directions for MIMO systems');
end

%   Default execution parameters
Def.maxiter = 50; 
Def.tol     = 1e-3; 
Def.type    = ''; %'stab'
Def.verbose = 0; % text output durint iteration?
Def.stopCrit= 'combAny'; %'s0', 'sysr', 'combAll', 'combAny'

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
r = length(s0);

s0old = s0;
s0 = cplxpair(s0);
if exist('Rt','var') && ~isempty(Rt)
    [~,cplxSorting] = ismember(s0,s0old); 
    Rt = Rt(:,cplxSorting);
else
    Rt = ones(sys.m,r);
end
clear s0old

% Initialize variables
sysr = sss([],[],[]);

if nargout > 7
    s0Traj = zeros(1,r,Opts.maxiter);
    RtTraj = zeros(sys.m,r,Opts.maxiter);

    s0Traj(:,:,1) = s0;
    RtTraj(:,:,1) = Rt;
end

%% Compute observability Gramian

% R = lyapchol(sys.A,sys.B,sys.E); %P = R'*R;
L = lyapchol(sys.A'/sys.E', sys.C'); Q = L'*L;

%% Run iteration

k=0;
while true
    k=k+1; sysr_old = sysr;
    
    if nargout > 7
        s0Traj(:,:,k) = s0;
        RtTraj(:,:,k) = Rt;
    end
    
    %--------   CORE  --------------
    %a) RK
    if sys.isSiso
        [~, V, ~,B_, Rsylv] = rk(sys, s0,Opts);
    else
        [~, V, ~, B_, Rsylv] = rk(sys, s0, Rt,Opts);
    end 
    %b) Left projection matrix
    W = Q*V/(V.'*Q*V);
    %c) ROM
    sysr = sss(W.'*sys.A*V, W.'*sys.B, sys.C*V, sys.D, W.'*sys.E*V);
    
    %d) Update of the reduction parameters
    s0_old=s0; if ~sys.isSiso, Rt_old = Rt; end
    if sys.isMimo
        [~, D, Y] = eig(sysr);
        Rt = full((Y.'*sysr.B).');
        s0 = full(-diag(D).');
        % make sure real shifts have real directions
        idx = find(imag(s0)==0); Rt(:,idx) = real(Rt(:,idx));
    else
        s0 = -eig(sysr)';
    end
    if strcmp(Opts.type,'stab')
        % mirror shifts with negative real part
        s0 = s0.*sign(real(s0));
    end   
    %-------------------------------
    
    [stop, stopCrit] = stoppingCriterion(s0,s0_old,sysr,sysr_old,Opts);
    if Opts.verbose
        fprintf('ISRK step %03u - Convergence: %s \n', ...
            k, sprintf('% 3.1e', stopCrit));
    end
    if stop || k>= Opts.maxiter
        s0 = s0_old; % function return value
        if ~sys.isSiso, Rt = Rt_old; end
        
        if nargout > 10
            % keep only what has been computed
            s0Traj = s0Traj(:,:,1:k); RtTraj = RtTraj(:,:,1:k);
        end
        break
    end      
end
if ~Opts.verbose %display at least the last value
    fprintf('ISRK step %03u - Convergence (%s): %s \n', ...
            k, Opts.stopCrit, sprintf('% 3.1e', stopCrit));
end
if k>=Opts.maxiter
    warning('ISRK:no_converged', ['ISRK has not converged after ' num2str(k) ' steps.']);
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

