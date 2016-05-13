function [sysr, V, W, s0] = ISRK(sys, s0, maxiter, epsilon, Qm) 
% ISRK Algorithm for Model Order Reduction by Krylov Subspace Methods
% ------------------------------------------------------------------
% [sysr, V, W, s0] = ISRK(sys, s0, maxiter, epsilon, Qm)
% Inputs:       * sys: an sss-object containing the LTI system
%               * s0: initial expansion points
%               * maxiter: maximum number of iterations
%               * epsilon: maximum step size for convergence
%               * Qm: Matrix definig direction of projection: W:=Qm*V
% Outputs:      * sysr: reduced system
%               * V: Projection matrix spanning a Krylov subspac
%               * W: Projection matrix
%               * s0: vector containing optimal expansion points
% ------------------------------------------------------------------
% For further explanation of the notation of s0, please refer to
% <a href="matlab:help RK">help RK</a>.
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de)
%               Thomas Wolf  (thomas.wolf@tum.de)
% Last Change:  20 Mar 2012
% ------------------------------------------------------------------

if ~exist('maxiter', 'var') || isempty(maxiter)
    maxiter=50;
end
if ~exist('epsilon', 'var') || isempty(epsilon)
    epsilon=1e-3;
elseif epsilon<=0 || ~isreal(epsilon)
    error('epsilon must be a real positive number.');
end

s0 = s0_vect(s0);

% Calculation of the Gramian if not known yet.
if ~exist('Qm', 'var')
    % Is Observability Gramian available?
    if isempty(sys.ObsGramChol)
        if isempty(sys.ObsGram)
           % No, it is not. Solve Lyapunov equation. 
           try
                if sys.is_dae
                    L = lyapchol(sys.A', sys.C', sys.E');
                else
                    L = lyapchol(sys.A', sys.C');
                end
                sys.ObsGramChol = sparse(L);
           catch ex
                warning(ex.message, 'Error in lyapchol. Trying without Cholesky factorization...')
                if sys.is_dae
                    sys.ObsGram = lyap(full(sys.A'), full(sys.C'*sys.C), [], full(sys.E'));
                else
                    sys.ObsGram = lyap(full(sys.A'), full(sys.C'*sys.C));
                end
%                 try
%                     L = chol(sys.ObsGram);
%                 catch ex
%                     myex = MException(ex2.identifier, ['System seems to be unstable. ' ex2.message]);
%                     throw(myex)
%                 end
           end

           % store system
           if inputname(1)
               assignin('caller', inputname(1), sys);
           end
        end
    end
end

tic

sysr=sss([],[],[]);

s0_traj = zeros(maxiter+2, length(s0));
s0_traj(1,:) = s0;

k=0;
while true
    k=k+1;
    
    sysr_old = sysr;
    s0_old   = s0;
    
    % Calculation of the Krylov subspace
    [~, V] = RK(sys, s0);
    
    % Perform reduction
    if exist('Qm', 'var')
        W = Qm*V;
    elseif ~isempty(sys.ObsGramChol)
        W = sys.ObsGramChol'*(sys.ObsGramChol*(sys.E*V));
    else
        W = sys.ObsGram*(sys.E*V);
    end
    sysr = sss(W'*sys.A*V, W'*sys.B, sys.C*V, sys.D, W'*sys.E*V);

    s0 = -eig(sysr)';
    s0(isnan(s0)) = 0;
    s0 = cplxpair(s0);
    
    
    % Analyse
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
        disp(['ISRK step ' num2str(k,'%03u') 9 ': Convergence '  num2str(es0, '%3.1e'), ', ' num2str(e, '%3.1e')]);
    end
    
    if es0 < epsilon || e < epsilon || k >= maxiter
        s0 = s0_old; % function return value
        s0_traj = s0_traj(1:(k+1),:);
        break
    end
end

assignin('base', 'ISRK_s0', s0_traj);

if sys.n<=500
    disp(['ISRK step ' num2str(k,'%03u') 9 ': Convergence '  num2str(es0, '%3.1e'), ', ' num2str(e, '%3.1e')]);
end

disp(['ISRK exit. Elapsed time: ' num2str(toc) 's']);

if k>=maxiter
    warning('ISRK:no_converged', ['ISRK has not converged after ' num2str(k) ' steps.']);
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

