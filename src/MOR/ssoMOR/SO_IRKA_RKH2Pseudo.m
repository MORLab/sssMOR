function [V,Mr,Kr,Dr,Br,iter] = SO_IRKA_RKH2Pseudo(tol,nr,nmax,M,K,D,B,M_ini,K_ini,D_ini,B_ini)
% SO-IRKA mit Residuenkorrektur nach dem H2-Pseudo-Optimum
%
% Input    -      tol:     Tolernz der Hauptschleife    
%                  nr:     Reduktionsordnung
%                nmax:     Maximale Anzahl der Iterationen
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
%                iter:     Anzahl der Iterationen

numfree = size(B,1);
n = 0;
[EP,Lb_n] = Neueptv(M_ini,D_ini,K_ini,B_ini);
V = Neuvr(EP,Lb_n);
es = 1;
while n <=nmax && es > tol
    M_n = V'*M*V;
    K_n = V'*K*V;
    D_n = V'*D*V;
    B_n = V'*B; 
    [EP_n,Lb_n] = Neueptv(M_n,D_n,K_n,B_n);
    V = Neuvr(EP_n,Lb_n);
    es = norm(EP_n-EP)/norm(EP); 
    EP = EP_n;
    n = n+1;
    if n == 3
        tol = tol*1e2;
    end
    if n == 5
        tol = tol*1e1;
    end
    if n == 7 && tol < 1e-8
        tol = tol*1e1;
    end
end

iter = n;
Mr = V'*M*V;
Kr = V'*K*V;
Dr = V'*D*V;
Br = V'*B;

    function [EP_n,Lb_n] = Neueptv(M_n,D_n,K_n,B_n)
        n0 = size(M_n,1);  
        g = sqrt(norm(K_n)/norm(M_n));
        d = 2/(norm(K_n)+g*norm(D_n));
        M_sc = g^2*d*M_n;
        D_sc = g*d*D_n;
        K_sc = d*K_n;
        E_2 = [eye(n0) zeros(n0);
            zeros(n0) M_n]; 
        E_2sc = [eye(n0) zeros(n0);
            zeros(n0) M_sc];
        A_2sc = [zeros(n0) eye(n0);
            -K_sc -D_sc];
        [EVr,EW,EVl] = eig(A_2sc,E_2sc,'qz');
        EW = g*diag(EW);
        [~,Index]=sort(real(EW),'descend');
        EW_sort = EW(Index,1);
        EVr_sort = EVr(:,Index);
        EVl_sort = EVl(:,Index);
        EVr_sort = [EVr_sort(1:n0,:);g*EVr_sort(n0+1:2*n0,:)];
        EVl_sort = [(1/g)*EVl_sort(1:n0,:);d*EVl_sort(n0+1:2*n0,:)];
        Lambda_E = EVl_sort'*E_2*EVr_sort;
        EP_n = - EW_sort(1:nr);
        lb = (Lambda_E\(EVl_sort'*([zeros(size(B_n));B_n]))).';
        Lb_n = lb*[eye(nr);zeros(2*n0-nr,nr)];
    end

    function V_n = Neuvr(EP_n,Lb_n)
        nnbc = size(Lb_n,1);
        nnr = nr/2;
        Ho = zeros(numfree,nnr*nnbc);
        Hrb = zeros(nnbc,nnr);
        Hrc = zeros(nnbc,nnr);
        V_n = zeros(numfree,nr);
        EPx = EP_n(1:2:nr);
        EPxx = [EP_n(1:2:nr-1);EP_n(2:2:nr)];
        Lb_n = Lb_n(:,1:2:nr-1);
        for i = 1:nnr
            Lb_n(:,i) = Lb_n(:,i)/norm(Lb_n(:,i));
        end
        for i = 1:nnr
            K_sigma = EPx(i)^2*M+EPx(i)*D+K;
            [L,U,P,Q,R] = lu(K_sigma);
            Ho(:,(i-1)*nnbc+1:i*nnbc) = Q*(U\(L\(P*(R\(B)))));
        end
        M11 = zeros(nnr,nnr);
        M21 = zeros(nnr,nnr);
        N11 = zeros(nnr,nnr);
        N21 = zeros(nnr,nnr);
        kmax = 50;
        e_lb = zeros(nnr,1);
        cm = 1;
        emax = 1;
        k = 1;
        while (k <= kmax) && (emax > 1e-4) && (cm < 1e4)
            for i = 1:nnr
                temp = B'*Ho(:,(i-1)*nnbc+1:i*nnbc);
                Hrb(:,i) = temp*Lb_n(:,i);
            end
            Hrb = conj(Hrb);
            for j = 1:nnr
                for i = 1:nnr
                    M11(i,j) = Lb_n(:,i).'*conj(Lb_n(:,j))/(conj(EPxx(j))+EPxx(i));
                    M21(i,j) = Lb_n(:,i)'*conj(Lb_n(:,j))/(EPxx(j)+EPxx(i));
                end
            end
            MM = [real(M11)+real(M21) imag(M11)+imag(M21);
                -imag(M11)+imag(M21) real(M11)-real(M21)];
            cm = cond(MM);
            Lctemp = ([real(Hrb) imag(Hrb)])/MM;
            Lcx = complex(Lctemp(:,1:nnr),Lctemp(:,nnr+1:2*nnr));
            for i = 1:nnr
                Lcx(:,i) = Lcx(:,i)/norm(Lcx(:,i));
            end
            for i = 1:nnr
                temp = Ho(:,(i-1)*nnbc+1:i*nnbc).'*B;
                Hrc(:,i) =  temp*Lcx(:,i);
            end
            Hrc = conj(Hrc);
            for j = 1:nnr
                for i = 1:nnr
                    N11(i,j) = Lcx(:,i).'*conj(Lcx(:,j))/(conj(EPxx(j))+EPxx(i));
                    N21(i,j) = Lcx(:,i)'*conj(Lcx(:,j))/(EPxx(j)+EPxx(i));
                end
            end
            NN = [real(N11)+real(N21) imag(N11)+imag(N21);
                -imag(N11)+imag(N21) real(N11)-real(N21)];
            Lbtemp = ([real(Hrc) imag(Hrc)])/NN;
            Lbx = complex(Lbtemp(:,1:nnr),Lbtemp(:,nnr+1:2*nnr));
            for i = 1:nnr
                Lbx(:,i) = Lbx(:,i)/norm(Lbx(:,i));
                e_lb(i) = Lbx(:,i)'*Lb_n(:,i)/norm(Lbx(:,i))/norm(Lb_n(:,i));
            end
            emax = 1-min(real(e_lb));
            Lb_n = Lbx;
            k = k+1;
        end
        for i = 1:nnr
            V_temp =  Ho(:,(i-1)*nnbc+1:i*nnbc)*Lb_n(:,i);
            V_n(:,2*i-1:2*i) = [real(V_temp) imag(V_temp)];
        end
        [V_n,~] = qr(V_n,0);
    end
end