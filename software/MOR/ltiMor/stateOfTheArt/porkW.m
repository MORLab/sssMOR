function [Ar,Br,Cr,Er] = porkW(W,S_W,Brt,B)
% porkW - Pseudo-Optimal Rational (Output) Krylov PORK [Wolf et al. 2013]
%
% Syntax: 
%       [Ar,Br,Cr,Er] = porkW(W,S_W,Brt,B)
%
% Description:
%
% Input Arguments:
%       *Required Input  Arguments:* 
%       -W,S_W,Brt:      solution of  W.'*A - S_W*W.'*E - Brt*C = 0
%       -B:              HFM input matrix
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

Pr_c = lyapchol(-S_W, Brt);
Cr = -Brt.'/Pr_c/Pr_c.';
Ar = S_W+Brt*Cr;
Br = W.'*B;
Er = eye(size(Ar));