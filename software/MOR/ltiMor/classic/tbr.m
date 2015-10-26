function [sysr, varargout] = tbr(sys, varargin)
% TBR - Performs model order reduction by the Truncated Balanced Realization
%
% Syntax:
%       [sysr, varargout] = TBR(sys, varargin)
%
%
% Inputs:
%       -sys:   an sss-object containing the LTI system
%       -q:     (opt.) order of reduced system
%
%
% Outputs:
%       -sysr:  reduced system
%       -V,W:   (opt.) projection matrices (only if q is given!)
%       -hsv:   Hankel singular values
%
%// Note: If no q is given, the balancing transformation and calculation of the
%// Hankel Singular Values is performed without subsequent model reduction.
%
%
% Examples:
%       No examples
%
%
%------------------------------------------------------------------
%   This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State Space, Model Order 
%   Reduction and System Analysis Toolbox developed at the Chair of 
%   Automatic Control, Technische Universitaet Muenchen. For updates 
%   and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
%   More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Sylvia Cremer, Rudy Eid
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  07 Feb 2011
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

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
            if sys.isdescriptor
                sys.ConGramChol = lyapchol(sys.A,sys.B,sys.E);
            else
                sys.ConGramChol = lyapchol(sys.A,sys.B);
            end
            R = sys.ConGramChol;
        catch ex
            warning(ex.message, 'Error in lyapchol. Trying without Cholesky factorization...')
            if sys.isdescriptor
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
            if sys.isdescriptor
                L = lyapchol(sys.A'/sys.E', sys.C');
            else
                L = lyapchol(sys.A',sys.C');
            end
            sys.ObsGramChol = sparse(L);
        catch ex
            warning(ex.message, 'Error in lyapchol. Trying without Cholesky factorization...')
            if sys.isdescriptor
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
