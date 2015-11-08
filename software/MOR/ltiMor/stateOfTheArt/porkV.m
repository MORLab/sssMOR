function [Ar,Br,Cr,Er] = porkV(V,S_V,Crt,C)
% porkV - Pseudo-Optimal Rational Krylov (Input)
%
% Syntax: 
%       [Ar,Br,Cr,Er] = porkV(V,S_V,Crt,C)
%
% Description:
%       TODO
% 
% Input Arguments:
%       *Required Input Arguments:*
%       -V,S_V,Crt:      solution of  A*V - E*V*S_V - B*Crt = 0
%       -C:              HFM output matrix
%
% Output Arguments: 
%       -Ar,Br,Cr,Er:    ROM matrices
% 
% Examples:
%       TODO
% 
% See Also: 
%       porkW, spark, rk
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

Qr_c = lyapchol(-S_V', Crt');
Br = -Qr_c\(Qr_c'\Crt');
Ar = S_V+Br*Crt;
Cr = C*V;
Er = eye(size(Ar));
