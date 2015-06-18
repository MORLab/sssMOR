function [J, g, H] = CostFunctionH2Bound(A, B, E, p)
% Cost Functional for H2 Error Bound Minimization, by A. Kohl
%   Input:  A,B,E: HFM matrices; 
%           p:       parameter vector [a,b]; 
%   Output: cost functional J; gradient g; Hessian H
% $\MatlabCopyright$

a = p(1); b = p(2); s0 = p(1)+[1 -1]*sqrt(p(1)^2-p(2));
As1 = A-s0(1)*E; As2 = A-s0(2)*E; L_E = chol(E);
r1 = As1\B; r2 = As1\(E*r1); r3 = As1\(E*r2);
LeB = L_E'\B;    
LeB_ = LeB + L_E*((r1 + (As2\(B+E*r1*2*a)))*2*a );
J = LeB_'*LeB_ - LeB'*LeB;
if (nargout==1), return; end

lAr1 = E*(As2\(A*r1)); lAr2 = E*(As2\(E*(As2\(A*r2))));
lAr3 = E*(As2\(E*(As2\(E*(As2\(A*r3))))));
lEr1 = E*(As2\(E*r1)); lEr2 = E*(As2\(E*(As2\(E*r2))));
lAEAr2 = E*(As2\(E*(As2\(A *(L_E\(L_E'\(A*r2)))))));
tmp = A *(L_E\(L_E'\(A*r3)));
lAEAr3 = E*(As2\(E*(As2\(E*(As2\tmp)))));
lAEAEAr3 = E*(As2\(E*(As2\(E*(As2\(A*(L_E\(L_E'\tmp))))))));
g_Ba = 4*a*lAEAr2 + 8*a^2*lAr2 - 4*a*b*lEr2 + 4*lAr1 + 4*a*lEr1;  g_Bb = -4*a*lAr2;

LeB_a = L_E'\g_Ba; LeB_b = L_E'\g_Bb; g = 2*real([LeB_'*LeB_a, LeB_'*LeB_b]);
g_Baa = 16*a*lAEAEAr3+32*a^2*lAEAr3-16*a*b*lAr3+12*lAEAr2+24*a*lAr2-4*b*lEr2+4*lEr1;
g_Bbb = 8*a*lAr3; g_Bab = -4*lAr2 -16*a*lAEAr3;

H = [LeB_a'*LeB_a + LeB_'*(L_E'\g_Baa), LeB_a'*LeB_b + LeB_'*(L_E'\g_Bab); ...
     0, LeB_b'*LeB_b + LeB_'*(L_E'\g_Bbb)];
H(2,1)=H(1,2); H=2*real(H);
end