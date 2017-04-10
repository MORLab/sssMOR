function [V,Mr,Kr,Dr,Br] = SO_TRKM(nr,M,K,D,B,M_ini,K_ini,D_ini,B_ini)
% Tangential-Rational-Krylov-Unterraum um Mehrpunkte f¨¹r
% Systeme zweiter Ordnung (SO-TRKM)
%
% Input    -       nr:     Reduktionsordnung
%                   M:     Massenmatrix des Originalsystems
%                   K:     Steifigkeitsmatrix des Originalsystems
%                   D:     D?mpfungsmatrix des Originalsystems
%                   B:     Eingangsmatrix des Originalsystems
%               M_ini:     Massenmatrix des initialen ROMs
%               K_ini:     Steifigkeitsmatrix des initialen ROMs
%               D_ini:     D?mpfungsmatrix des initialen ROMs
%               B_ini:     Eingangsmatrix des initialen ROMs
% 
% Output    -       V:     Rechte Projektionsmatrix
%                  Mr:     Massenmatrix des ROMs
%                  Kr:     Steifigkeitsmatrix des ROMs
%                  Dr:     D?mpfungsmatrix des ROMs
%                  Br:     Eingangsmatrix des ROMs

numfree = size(B,1);
[EP,Lb] = Neueptv(M_ini,D_ini,K_ini,B_ini);
V = zeros(numfree,nr);
i = 1;
while i <= nr
    K_sigma = EP(i)^2*M+EP(i)*D+K;
    [L,U,P,Q,R] = lu(K_sigma);
    V_temp = Q*(U\(L\(P*(R\(B*Lb(:,i))))));
    if ~isreal(EP(i))
        V(:,i:i+1) = [real(V_temp) imag(V_temp)];
        i = i+2;
    else
        V(:,i) = V_temp;
        i = i+1;
    end
    [V(:,1:i-1),~] = qr(V(:,1:i-1),0);
end      
Mr = V'*M*V;
Kr = V'*K*V;
Dr = V'*D*V;
Br = V'*B;

    function [EPn,Lbn] = Neueptv(M_ini,D_ini,K_ini,B_ini)
        n0 = size(M_ini,1);
        g = sqrt(norm(K_ini)/norm(M_ini));
        d = 2/(norm(K_ini)+g*norm(D_ini));
        M_sk = g^2*d*M_ini;
        D_sk = g*d*D_ini;
        K_sk = d*K_ini;
        E_2 = [eye(n0) zeros(n0);
            zeros(n0) M_ini];
        E_2sc = [eye(n0) zeros(n0);
            zeros(n0) M_sk];
        A_2sc = [zeros(n0) eye(n0);
            -K_sk -D_sk];
        [EVr,EW,EVl] = eig(A_2sc,E_2sc,'qz');
        EW = g*diag(EW);
        [~,Index]=sort(real(EW),'descend');
        EW_sort = EW(Index,1);
        EVr_sort = EVr(:,Index);
        EVl_sort = EVl(:,Index);
        EVr_sort = [EVr_sort(1:n0,:);g*EVr_sort(n0+1:2*n0,:)];
        EVl_sort = [(1/g)*EVl_sort(1:n0,:);d*EVl_sort(n0+1:2*n0,:)];
        Lamda_E = EVl_sort'*E_2*EVr_sort;
        EPn = - EW_sort(1:nr);
        lb = (Lamda_E\(EVl_sort'*([zeros(size(B_ini));B_ini]))).';
        Lbn = lb*[eye(nr);zeros(2*n0-nr,nr)];
    end
end