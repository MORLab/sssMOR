function [sysr, varargout] = tbr(sys, varargin)
% Performs model order reduction by the Truncated Balanced Realization
% ------------------------------------------------------------------
% [sysr, varargout] = TBR(sys, varargin)
% Inputs:       * sys: an sss-object containing the LTI system
%    [optional] * q: order of reduced system
% Output:       * sysr: reduced system
%    [optional] * V, W: projection matrices (only if q is given!)
%               * hsv: Hankel singular values
% ------------------------------------------------------------------
% If no q is given, the balancing transformation and calculation of the
% Hankel Singular Values is performed without subsequent model reduction.
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Sylvia Cremer, Rudy Eid
% Last Change:  07 Feb 2011
% ------------------------------------------------------------------

% if sys.is_dae
% % if E-speye(size(E)) == sparse(size(E,1), size(E,1))
%     %mit E^-1 durchmultipliziert
% %     S = lyapchol(A,B,E);
% %     R = lyapchol(transpose(E\A),C');
%     %Zustandstrafo
%     if ~isfield(sys,'P_chol')
%         S = lyapchol(A,B,E);
% %         S = lyapchol(A/E,B); % P=S'*S
%         sys.P_chol=sparse(S);
%     else
%         S=full(sys.P_chol);
%     end
%     if ~isfield(sys,'Q_chol')
%         R = lyapchol(transpose(E\A),C');
% %         R = lyapchol(A',C',E'); % Q=R'*R
%         sys.Q_chol=sparse(R);
%     else
%         R=full(sys.Q_chol);
%     end
% else
%     if ~isfield(sys,'P_chol')
%         S = lyapchol(A,B); % P=S'*S
%         sys.P_chol=sparse(S);
%     else
%         S=full(sys.P_chol);
%     end
%     if ~isfield(sys,'Q_chol')
%         R = lyapchol(A',C'); % Q=R'*R
%         sys.Q_chol=sparse(R);
%     else
%         R=full(sys.Q_chol);
%     end 
% end
    
% Is Controllability Gramian available?
if isempty(sys.ConGramChol)
    if isempty(sys.ConGram)
        % No, it is not. Solve Lyapunov equation.
        try
            if sys.isDescriptor
                sys.ConGramChol = lyapchol(sys.A,sys.B,sys.E);
            else
                sys.ConGramChol = lyapchol(sys.A,sys.B);
            end
            R = sys.ConGramChol;
        catch ex
            warning(ex.message, 'Error in lyapchol. Trying without Cholesky factorization...')
            if sys.isDescriptor
                try
                    sys.ConGram = lyap(sys.A, sys.B*sys.B', [], sys.E);
                catch ex2
                    warning(ex2.message, 'Error solving Lyapunov equation. Premultiplying by E^(-1)...')
                    tmp = sys.E\sys.B;
                    sys.ConGram = lyap(sys.E\sys.A, tmp*tmp');
                end
            else
                sys.ConGram = lyap(full(sys.A), full(sys.B*sys.B'));
            end
            try
                R = chol(sys.ConGram);
            catch ex2
                myex = MException(ex2.identifier, ['System seems to be unstable. ' ex2.message]);
                throw(myex)
            end
        end
    else
        R = chol(sys.ConGram);
    end
else
    R = sys.ConGramChol;
end


% Is Observability Gramian available?
if isempty(sys.ObsGramChol)
    if isempty(sys.ObsGram)
        % No, it is not. Solve Lyapunov equation. 
       try
            if sys.isDescriptor
                L = lyapchol(sys.A'/sys.E', sys.C');
            else
                L = lyapchol(sys.A',sys.C');
            end
            sys.ObsGramChol = sparse(L);
        catch ex
            warning(ex.message, 'Error in lyapchol. Trying without Cholesky factorization...')
            if sys.isDescriptor
                sys.ObsGram = lyap(sys.A'/sys.E', sys.C'*sys.C);
            else
                sys.ObsGram = lyap(full(sys.A'), full(sys.C'*sys.C));
            end
            try
                L = chol(sys.ObsGram);
            catch ex
                myex = MException(ex2.identifier, ['System seems to be unstable. ' ex2.message]);
                throw(myex)
            end
       end
    else
        L = chol(sys.ObsGram);
    end
else
    L = sys.ObsGramChol;
end

% calculate balancing transformation and Hankel Singular Values
% M = L*R';
% [K,S] = svd(M*M');
% hsv = sqrt(diag(S));
% sys.HankelSingularValues = real(hsv);
% sys.T_bal = diag(sqrt(hsv))\K'*L/sys.E;
% sys.T_bal_inv = L\K*diag(sqrt(hsv));

[K,S,M]=svd(R*L');
hsv = diag(S);
sys.HankelSingularValues = real(hsv);
sys.T_bal_inv = R'*K/diag(sqrt(hsv));
sys.T_bal = diag(sqrt(hsv))\M'*L/sys.E;


% store system
if inputname(1)
    assignin('caller', inputname(1), sys);
end

if nargin==1
    return
end

% if MOR is to be performed, calculate V, W and reduced system
q=varargin{1};

V = sys.T_bal_inv(:,1:q);
W = sys.T_bal(1:q,:)';

sysr = sss(W'*sys.A*V, W'*sys.B, sys.C*V, sys.D, W'*sys.E*V);
varargout{1} = V;
varargout{2} = W;
varargout{3} = real(hsv);
