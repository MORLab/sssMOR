function [mu, L_S, P_S] = numericalAbscissa(A,E)
% numericalAbscissa - Compute generalized numerical abscissa
%
% Syntax:
%       [mu, L_S, P_S] = numericalAbscissa(A,E)
%
% Description:
%	   TODO
%
% Input Arguments:
%	   A, E:     HFM matrices
%
% Output Arguments:
%	   mu:       generalized numerical abscissa mu_E(A)
%       L_S, P_S: Cholesky decomposition of S=-A-A'
%
% See Also: 
%       TODO
%
% References:
%	   TODO
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
% Authors:      TODO
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  03 Mai 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

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