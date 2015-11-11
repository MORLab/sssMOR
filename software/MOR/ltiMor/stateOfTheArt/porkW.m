function varargout = porkW(W,S_W,Brt,B)
% PORKW - Pseudo-Optimal Rational Krylov (Output)
%
% Syntax: 
%       [Ar,Br,Cr,Er] = porkW(W,S_W,Brt,B)
%
% Description:
%       This function implements the pseudo-optimal rational Krylov
%       algorithm introduced by Wolf and Panzer [1,2].
%
%       Given a projection matrix W spanning an output Krylov subspace and
%       the respective matrices from the Sylvester equation
%
%          $ A^T W - E^T W S_W^T - C^T Brt^T = 0$
%
%       this function computes the reduced order matrices corresponding to
%       the H2-pseudo-optimal reduced order model, i.e. a model
%       interpolating the original according to (W,S_W,Brt) and having
%       eigenvalues as mirror images of the shifts.
%
% Input Arguments:
%       *Required Input  Arguments:* 
%       -W,S_W,Brt:      solution of  W.'*A - S_W*W.'*E - Brt*C = 0
%       -B:              input matrix of the original model
%
% Output Arguments: 
%       - sysrPO:         Pseudo-optimal reduced order model 
%       - Ar,Br,Cr,Er:    ROM matrices
%
% Examples:
%       Following code computes an H2-pseudo-optimal reduced order model
%       with an output Krylov subspace
%> sys = loadSss('build');
%> s0 = -eigs(sys,4,'sm').';
%> [sysr, ~, W] = rk(sys,[],s0);
%> [Brt, S] = getSylvester(sys, sysr, W, 'W');
%> sysrPO = porkW(V,S,Brt',sys.B)
% 
% See Also: 
%       porkV, spark, rk, getSylvester
%
% References:
%       * *[1] Wolf (2014)*, H2 Pseudo-Optimal Moder Order Reduction
%       * *[2] Panzer (2014)*, Model Order Reduction by Krylov Subspace Methods
%              with Global Error Bounds and Automatic Choice of Parameters
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
% Authors:      Thomas Wolf, Heiko Panzer
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Computation
Pr_c = lyapchol(-S_W, Brt);
Cr = -Brt.'/Pr_c/Pr_c.';
Ar = S_W+Brt*Cr;
Br = W.'*B;
Er = eye(size(Ar));

%% Preparing output
if nargout == 1
    varargout{1} = sss(Ar,Br,Cr,zeros(size(Cr,1),size(Br,2)),Er);
else
    varargout{1} = Ar;
    varargout{2} = Br;
    varargout{3} = Cr;
    varargout{4} = Er;
end
