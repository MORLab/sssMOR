function [sysr, V, W, s0, Rt, Lt, B_, Sv, Rv, C_, Sw, Lw, kIter, s0Traj, RtTraj, LtTraj] = irka(sys, varargin) 
% IRKA - Iterative Rational Krylov Algorithm
%
% Syntax:
%       sysr                            = IRKA(sys, q)
%       sysr                            = IRKA(sys, s0)
%       sysr                            = IRKA(sys, s0, Rt, Lt)
%       sysr                            = IRKA(sys, s0,..., Opts)
%       [sysr, V, W]                    = IRKA(sys, s0,... )
%       [sysr, V, W, s0]                = IRKA(sys, s0,... )
%       [sysr, V, W, s0, Rt, Lt]        = IRKA(sys, s0,... )
%       [sysr, V, W, s0, Rt, Lt, B_, Sv, Rv, C_, Sw, Lw, kIter] = IRKA(sys, s0,... )
%       [sysr, V, W, s0, Rt, Lt, B_, Sv, Rv, C_, Sw, Lw, kIter,  s0Traj, RtTraj, LtTraj] = IRKA(sys, s0,... )
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
%       The algorithm can be called either with a reduced order q, in which
%       case all initial shifts are set at the origin (and tangential
%       directions to ones), or with a prescribed initial set of shifts s0
%       (and tangential directions Rt,Lt).
%
% Input Arguments:  
%       *Required Input Arguments:*
%       -sys:			full order model (sss)
%
%       *Optional Input Arguments:*
%       -q:             reduced order
%       -s0:			vector of initial shifts
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
%           -.(refer to *arnoldi* or *rk* for other options)
%
% Output Arguments:      
%       -sysr:              reduced order model (sss)
%       -V,W:               resulting projection matrices
%       -s0:                final choice of shifts
%       -Rt,Lt:             final choice of tangential directions
%       -B_,Sv,Rv:          matrices of the input Sylvester equation
%       -C_,Sw,Lw:          matrices of the output Sylvester equation
%       -kIter:             number of iterations
%		-s0Traj,RtTraj,..:  trajectory of all shifst and tangential directions for all iterations
%
% Examples:
%       This code computes an H2-optimal approximation of order 8 to
%       the benchmark model 'building'. One can use the function isH2opt to
%       verify if the necessary conditions for optimality are satisfied.
%
%> sys = loadSss('building')
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
% Last Change:  20 Jan 2017
% Copyright (c) 2016,2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Parse input and load default parameters
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin = varargin(1:end-1);
else
    Opts = struct();
end
if isempty(varargin)
        error('sssMOR:irka:inputarguments','Not enough input arguments');
elseif length(varargin)==1
    if length(varargin{1}) == 1
        % Reduced order specified
        s0 = zeros(1,varargin{1});
    else
        % Shifts specified
        s0 = varargin{1};
    end
    Rt = ones(sys.m,length(s0));
    Lt = ones(sys.p,length(s0));
else %MIMO
        %usage irka(sys,s0,Rt,Lt)
        s0 = varargin{1};
        Rt = varargin{2};
        Lt = varargin{3};
end
    
%% Parse the inputs
%   Default execution parameters
Def.maxiter = 50; 
Def.tol     = 1e-3; 
Def.type    = ''; %'stab', 'newton', 'restarted'
Def.verbose = 0; % text output durint iteration?
Def.stopCrit= 'combAny'; %'s0', 'sysr', 'combAll', 'combAny'
Def.suppressverbose = 0;

% create the options structure
if ~exist('Opts','var') || isempty(fieldnames(Opts))
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

% sort expansion points & tangential directions
s0old = s0;
s0 = cplxpair(s0);
if exist('Rt','var') && ~isempty(Rt)
    [~,cplxSorting] = ismember(s0,s0old); 
    Rt = Rt(:,cplxSorting);
    Lt = Lt(:,cplxSorting);
else
    Rt = ones(sys.m,r);
    Lt = ones(sys.p,r);
end
clear s0old

% Initialize variables
sysr = sss([],[],[]);

s0Traj = zeros(1,r,Opts.maxiter);
RtTraj = zeros(sys.m,r,Opts.maxiter);
LtTraj = zeros(sys.p,r,Opts.maxiter);

s0Traj(:,:,1) = s0;
RtTraj(:,:,1) = Rt;
LtTraj(:,:,1) = Lt;

%% IRKA iteration
kIter=0;
while true
    kIter=kIter+1; sysr_old = sysr;
    %   Reduction
    if sys.isSiso
        [sysr, V, W, B_, Sv, Rv, C_, Sw, Lw] = rk(sys, s0, s0,Opts);
    else
        [sysr, V, W, B_, Sv, Rv, C_, Sw, Lw] = rk(sys, s0, s0, Rt, Lt,Opts);
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
            warning('sssMOR:irka:tangDirIncompatible',...
                'Tangential directions corresponding to real shifts are complex')              
        end
        Rt(:,idx) = real(Rt(:,idx)); Lt(:,idx) = real(Lt(:,idx));
    else
        s0 = -eig(sysr)';
    end
    if strcmp(Opts.type,'stab')
        % mirror shifts with negative real part
        s0 = s0.*sign(real(s0));
    end

    s0Traj(:,:,kIter+1) = s0;
    RtTraj(:,:,kIter+1) = Rt; 
    LtTraj(:,:,kIter+1) = Lt;
    
    [stop, stopCrit] = stoppingCriterion(s0,s0_old,sysr,sysr_old,Opts);
    if Opts.verbose
        fprintf('IRKA step %03u - Convergence: %s \n', ...
            kIter, sprintf('% 3.1e', stopCrit));
    end
    if stop || kIter>= Opts.maxiter
        s0 = s0_old; % function return value
        if ~sys.isSiso, Rt = Rt_old; Lt = Lt_old; end
        % keep only what has been computed
        s0Traj = s0Traj(:,:,1:kIter);
        RtTraj = RtTraj(:,:,1:kIter); LtTraj = LtTraj(:,:,1:kIter);
        break
    end      
end

%%  Storing additional parameters
%Stroring additional information about the reduction in the object 
%containing the reduced model:
%   1. Define a new field for the Opts struct and write the information
%      that should be stored to this field
%   2. Adapt the method "parseParamsStruct" of the class "ssRed" in such a
%      way that the new defined field passes the check
Opts.originalOrder = sys.n;
if ~isfield(Opts,'orth'),   Opts.orth = '2mgs'; end
if ~isfield(Opts,'reorth'), Opts.reorth = 0; end
if ~isfield(Opts,'lse'),    Opts.lse = 'sparse'; end
if ~isfield(Opts,'dgksTol'),Opts.dgksTol = 1e-12; end
if ~isfield(Opts,'krylov'), Opts.krylov = 0; end
Opts.Rt = Rt;
Opts.Lt = Lt;
Opts.kIter = kIter;
Opts.s0Traj = s0Traj;
Opts.RtTraj = RtTraj;
Opts.LtTraj = LtTraj;
Opts.s0 = s0;

% Convert the reduced system to a ssRed-object
sysr = ssRed(sysr.A,sysr.B,sysr.C,sysr.D,sysr.E,'irka',Opts,sys);

%%  Finish execution
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
        if any(abs(s0))<1e-3
             stopCrit = norm((s0-s0_old), 1)/sysr.n;
        else
            stopCrit = norm((s0-s0_old)./s0, 1)/sysr.n;
        end      
        stop = stopCrit <= Opts.tol;
    case 'sysr' %reduced model convergence
        stopCrit = inf; %initialize in case the reduced model is unstable
        if all(real(eig(sysr))<0) && all(real(eig(sysr_old))<0)
                stopCrit=norm(sysr-sysr_old)/norm(sysr);
        end
        stop = stopCrit <= Opts.tol;
    case 'combAll'
        if any(abs(s0))<1e-3
             stopCrit = norm((s0-s0_old), 1)/sysr.n;
        else
            stopCrit = norm((s0-s0_old)./s0, 1)/sysr.n;
        end 
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
