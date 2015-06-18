function [Ar,Br,Cr,Er] = PORK_V(V,S_V,Crt,C)
% Pseudo-Optimal Rational (Input) Krylov PORK [Wolf et al. 2013]
%   Input:  V,S_V,Crt:      solution of  A*V - E*V*S_V - B*Crt = 0
%           C:              HFM output matrix
%   Output: Ar,Br,Cr,Er:    ROM matrices
% $\MatlabCopyright$

Qr_c = lyapchol(-S_V', Crt');
Br = -Qr_c\(Qr_c'\Crt');
Ar = S_V+Br*Crt;
Cr = C*V;
Er = eye(size(Ar));
