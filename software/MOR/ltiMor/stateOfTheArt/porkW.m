function [Ar,Br,Cr,Er] = PORK_W(W,S_W,Brt,B)
% Pseudo-Optimal Rational (Output) Krylov PORK [Wolf et al. 2013]
%   Input:  W,S_W,Brt:      solution of  W.'*A - S_W*W.'*E - Brt*C = 0
%           B:              HFM input matrix
%   Output: Ar,Br,Cr,Er:    ROM matrices
% $\MatlabCopyright$

Pr_c = lyapchol(-S_W, Brt);
Cr = -Brt.'/Pr_c/Pr_c.';
Ar = S_W+Brt*Cr;
Br = W.'*B;
Er = eye(size(Ar));