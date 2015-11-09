function [Ar,Br,Cr,Er] = porkW(W,S_W,Brt,B)
% porkW - Pseudo-Optimal Rational Krylov (Output)
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
%           A.' W - E.' W S_W.' - C.' Brt.' = 0
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
%       -Ar,Br,Cr,Er:    ROM matrices
%
% Examples:
%       TODO
% 
% See Also: 
%       porkV, spark, rk
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

Pr_c = lyapchol(-S_W, Brt);
Cr = -Brt.'/Pr_c/Pr_c.';
Ar = S_W+Brt*Cr;
Br = W.'*B;
Er = eye(size(Ar));
