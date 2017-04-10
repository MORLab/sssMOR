function [V,Mr,Kr,Dr,Br] = CraigBampton(nr,M,K,D,B,interndof,nbc)
% Craig-Bampton-Verfahren
%
% Input    -       nr:     Reduktionsordnung
%                   M:     Massenmatrix des Originalsystems
%                   K:     Steifigkeitsmatrix des Originalsystems
%                   D:     D?mpfungsmatrix des Originalsystems
%                   B:     Eingangsmatrix des Originalsystems
%           interndof:     Innere FHG
%                 nbc:     Schnittstellen-FHG
% 
% Output    -       V:     Rechte Projektionsmatrix
%                  Mr:     Massenmatrix des ROMs
%                  Kr:     Steifigkeitsmatrix des ROMs
%                  Dr:     D?mpfungsmatrix des ROMs
%                  Br:     Eingangsmatrix des ROMs

numnbc = length(nbc);
num_schwinmod = nr - numnbc;
schnittdof = nbc;
B = [B(nbc,:);B(interndof,:)];
M_ii = M(interndof,interndof);
K_ii = K(interndof,interndof);
K_ib = K(interndof,schnittdof);
[Vd,~] = eigs(K_ii,M_ii,num_schwinmod, 'sm');
Vd = [zeros(numnbc, num_schwinmod);Vd];
[R,~,S] = chol(K_ii);
X = S*(R\(R'\(S'*K_ib)));
Vs = [eye(numnbc);-X];
V = [Vs Vd];
M = [M(schnittdof,schnittdof) M(schnittdof,interndof);
      M(interndof,schnittdof) M_ii];
K = [K(schnittdof,schnittdof) K(schnittdof,interndof);
      K_ib K_ii];
D = [D(schnittdof,schnittdof) D(schnittdof,interndof);
      D(interndof,schnittdof) D(interndof,interndof)];
Mr = V'*M*V;
Kr = V'*K*V;
Dr = V'*D*V;
Br = V'*B;
end