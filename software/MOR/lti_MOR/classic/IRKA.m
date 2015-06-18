function [sysr, V, W, s0, s0_traj] = IRKA(sys, s0, maxiter, epsilon,opts) 
% IRKA Algorithm for Model Order Reduction by Krylov Subspace Methods
% ------------------------------------------------------------------
% [sysr, V, W, s0] = IRKA(sys, [s0;q], maxiter, epsilon)
% Inputs:       * sys: an sss-object containing the LTI system
%               * s0: initial expansion point
%               * maxiter: maximum number of iterations
%               * epsilon: maximum step size for convergence
%               * opts: choose between original IRKA ('orig') and a 
%                       stability preserving version ('stab')
% Outputs:      * sysr: reduced system
%               * V, W: Projection matrices spanning Krylov subspaces
%               * s0: vector containing optimal expansion points
%               * s0_traj: complete history of s0 over iterations
% ------------------------------------------------------------------
% For further explanation of the notation of s0, please refer to
% <a href="matlab:help RK">help RK</a>.
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto (a.castagnotto@tum.de),
%               Heiko Panzer (heiko@mytum.de)
% Last Change:  24 Feb 2015
% ------------------------------------------------------------------

if ~exist('maxiter', 'var') || isempty(maxiter)
    maxiter=50;
end
if ~exist('epsilon', 'var') || isempty(epsilon)
    epsilon=1e-3;
elseif epsilon<=0 || ~isreal(epsilon)
    error('epsilon must be a real positive number.');
end
if ~exist('opts', 'var') || isempty(opts)
    opts = 'stab';
end

s0 = s0_vect(s0);
tic

sysr=sss([],[],[]);

s0_traj = zeros(maxiter+2, length(s0));
s0_traj(1,:) = s0;

k=0;
while true
    k=k+1;
    
    sysr_old = sysr;
    [sysr, V, W] = RK(sys, s0, s0); %@(x,y) (x'*y)

    s0_old=s0;
    s0 = -eig(sysr)';

    s0(isnan(s0)) = 0;
    s0 = cplxpair(s0);

    if strcmp(opts,'stab')
        % mirror shifts with negative real part
        s0 = s0.*sign(real(s0));
    end

    s0_traj(k+1,:) = s0;
%     plot([s0_old;s0], '-*'); hold on
    
    es0 = norm((s0-s0_old)./s0, 1)/sysr.n;

    e=inf;
    if all(real(s0)>0)
        warning off
        if isreal(norm(sysr))
            e=norm(sysr-sysr_old)/norm(sysr);
        end
        warning on
    end
    
    if sys.n>500
        disp(['IRKA step ' num2str(k,'%03u') 9 ': Convergence '  num2str(es0, '%3.1e'), ', ' num2str(e, '%3.1e')]);
    end
    
    if es0 < epsilon || k >= maxiter %  || e < epsilon
        s0 = s0_old; % function return value
        s0_traj = s0_traj(1:(k+1),:);
        break
    end
end

% assignin('base', 'IRKA_s0', s0_traj);

if sys.n<=500
    disp(['IRKA step ' num2str(k,'%03u') 9 ': Convergence '  num2str(es0, '%3.1e'), ', ' num2str(e, '%3.1e')]);
end

disp(['IRKA exit. Elapsed time: ' num2str(toc) 's']);

if k>=maxiter
    warning('IRKA:no_converged', ['IRKA has not converged after ' num2str(k) ' steps.']);
    return
end

end

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
end

