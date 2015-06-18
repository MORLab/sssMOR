function [mu, L_S, P_S] = NumericalAbscissa(A,E)
% Compute generalized numerical abscissa
%   Input:  A, E:     HFM matrices
%   Output: mu:       generalized spectral abscissa mu_E(A)
%           L_S, P_S: Cholesky decomposition of S=-A-A'
% $\MatlabCopyright$

p    = 20;		% number of Lanczos vectors
tol  = 1e-10;	% convergence tolerance
opts = struct('issym',true, 'p',p, 'tol',tol, 'v0',sum(E,2));

[~,e_E,~] = chol(sparse(E)); 
if e_E, error('E is not symmetric positiv definite.'); end

[L_S,e,P_S] = chol(sparse(-A-A')); % L_S'*L_S = P_S'*(-A-A')*P_S
if e, error('A is not strictly dissipative.'); 
else
    symA = @(x) P_S*(L_S\(L_S'\(P_S'*x*2)));
    mu   = -eigs(symA, size(A,1), E, 1, 0, opts);
end