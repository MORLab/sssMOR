function [Ar,Br,Cr,Er] = porkV(V,S_V,Crt,C)
% porkV - Pseudo-Optimal Rational (Input) Krylov PORK [Wolf et al. 2013]
%
% Syntax: 
%       [Ar,Br,Cr,Er] = porkV(V,S_V,Crt,C)
%
% Description:
% 
% Input Arguments:
%       *Required Input Arguments:*
%       -V,S_V,Crt:      solution of  A*V - E*V*S_V - B*Crt = 0
%       -C:              HFM output matrix
%
% Output Arguments: 
%       -Ar,Br,Cr,Er:    ROM matrices
% 
% References:
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
% Authors:      
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
