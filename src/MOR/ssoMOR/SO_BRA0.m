function [V,Mr,Kr,Dr,Br] = SO_BRA0(nr,M,K,D,B)
% Block-Rational-Arnoldi um Null f¨¹r Systeme zweiter
% Ordnung (SO-BRAE)
%
% Input    -       nr:     Reduktionsordnung
%                   M:     Massenmatrix des Originalsystems
%                   K:     Steifigkeitsmatrix des Originalsystems
%                   D:     D?mpfungsmatrix des Originalsystems
%                   B:     Eingangsmatrix des Originalsystems
% 
% Output    -       V:     Rechte Projektionsmatrix
%                  Mr:     Massenmatrix des ROMs
%                  Kr:     Steifigkeitsmatrix des ROMs
%                  Dr:     D?mpfungsmatrix des ROMs
%                  Br:     Eingangsmatrix des ROMs

numfree = size(B,1);
numnbc = size(B,2);
V = zeros(numfree,nr);
[R,~,S] = chol(K);
V0 = S*(R\(R'\(S'*B)));
[Q0,~] = qr(V0,0);
V(:,1:numnbc) = Q0;
i = 1;

while i+numnbc<=nr
    Vi = M*V(:,i);
    r = S*(R\(R'\(S'*Vi)));
    for j = 1:i+numnbc-1
        h = V(:,j)'*r;
        r = r-h*V(:,j);
    end
    h_neu = norm(r);
    V(:,i+numnbc) = r/h_neu;
    i = i+1;
end
Mr = V'*M*V;
Kr = V'*K*V;
Dr = V'*D*V;
Br = V'*B;
end