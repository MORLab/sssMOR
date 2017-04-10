function [V,Mr,Kr,Dr,Br,iter] = SO_IRKA_RKQuasiNewton(tol,nr,nmax,M,K,D,B,M_ini,K_ini,D_ini,B_ini)
% SO-IRKA mit Residuenkorrektur mit Quasi-Newton-Verfahren
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
[EP,Lb_n,Lc_n] = Neueptv(M_ini,D_ini,K_ini,B_ini);
V = Neuvr(EP,Lb_n,Lc_n);
es = 1;
while n <=nmax && es > tol
    M_n = V'*M*V;
    K_n = V'*K*V;
    D_n = V'*D*V;
    B_n = V'*B; 
    [EP_n,Lb_n,Lc_n] = Neueptv(M_n,D_n,K_n,B_n);
    V = Neuvr(EP_n,Lb_n,Lc_n);
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

    function [EP_n,Lb_n,Lc_n] = Neueptv(M_n,D_n,K_n,B_n)
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
        Lb_n = ([eye(nr) zeros(nr,2*n0-nr)]*(Lambda_E\(EVl_sort'*([zeros(size(B_n));B_n])))).';
        Lc_n = ([B_n' zeros(size(B_n'))])*EVr_sort*[eye(nr);zeros(2*n0-nr,nr)];
    end

    function V_n = Neuvr(EP,Lb_n,Lc_n)
        nnbc = size(Lb_n,1);
        nnr = nr/2;
        Ho = zeros(numfree,nnr*nnbc);
        V_n = zeros(numfree,nr);
        EPx = EP(1:2:nr);
        Lb_n = Lb_n(:,1:2:nr-1);
        Lc_n = Lc_n(:,1:2:nr-1);
        for i = 1:nnr
            K_sigma = EPx(i)^2*M+EPx(i)*D+K;
            [L,U,P,Q,R] = lu(K_sigma);
            Ho(:,(i-1)*nnbc+1:i*nnbc) = Q*(U\(L\(P*(R\(B)))));
        end
        Hobc = B'*Ho;
        z0 = zeros(2*nnr*nnbc,1);
        for i = 1:nnr
            z0(2*(i-1)*nnbc+1:2*i*nnbc) = [Lb_n(:,i);Lc_n(:,i)];
        end
        [~,gz0] = fobject(z0,Hobc,EPx);
        p0 = -gz0;
        alpha0 = linesearch(@fobject,z0,p0,Hobc,EPx);
        s0 = alpha0*p0;
        [~,gz1] = fobject(z0+s0,Hobc,EPx);
        y0 = gz1-gz0;
        H0 = real(y0'*s0)/(y0'*y0)*eye(2*nnr*nnbc);
        Hk = H0;
        gzk = gz0;
        zk = z0;
        k = 1;
        kmax = 100;
        while (norm(gzk) > 1e-5) && (k <= kmax)
            pk = -Hk*gzk;
            alphak = linesearch(@fobject,zk,pk,Hobc,EPx);
            sk = alphak*pk;
            [~,gzkn] = fobject(zk+sk,Hobc,EPx);
            yk = gzkn - gzk;
            rhok = 1/real(yk'*sk);
            Hk = ((eye(2*nnr*nnbc)-rhok*(sk*yk'))*Hk)*(eye(2*nnr*nnbc)-rhok*(yk*sk')) + rhok*(sk*sk');
            gzk = gzkn;
            zk = zk + sk;
            k = k + 1;
        end
        for i = 1:nnr
            V_temp =  Ho(:,(i-1)*nnbc+1:i*nnbc)*zk(2*(i-1)*nnbc+1:(2*i-1)*nnbc);
            V_n(:,2*i-1:2*i) = [real(V_temp) imag(V_temp)];
        end
        [V_n,~] = qr(V_n,0);
    end    
    
    function [fz,gzs] = fobject(z,Hobc,EPx)        
        nnr = size(EPx,1);
        nnbc = size(Hobc,1);
        gzs = zeros(2*nnr*nnbc,1);
        fz = 0;
        for i = 1:nnr
            temp_Ho = conj(Hobc(:,(i-1)*nnbc+1:i*nnbc));
            temp_Jbcx = zeros(nnbc,nnbc);
            for k = 1:nnr
                temp_phi = z((2*k-1)*nnbc+1:2*k*nnbc)*z(2*(k-1)*nnbc+1:(2*k-1)*nnbc).';
                temp_Jbcx = temp_Jbcx + conj(temp_phi/(EPx(i)+EPx(k)))+temp_phi/(conj(EPx(i))+EPx(k));
            end
            fz = fz + 2*real((z((2*i-1)*nnbc+1:2*i*nnbc)'*(-2*temp_Ho+temp_Jbcx))*conj(z(2*(i-1)*nnbc+1:(2*i-1)*nnbc)));
            Jbx = 2*z((2*i-1)*nnbc+1:2*i*nnbc)'*(-temp_Ho+temp_Jbcx);
            Jcx = 2*(-temp_Ho+temp_Jbcx)*conj(z(2*(i-1)*nnbc+1:(2*i-1)*nnbc));
            gzs(2*(i-1)*nnbc+1:2*i*nnbc) = [Jbx.';Jcx];
        end
    end

    function alphas = linesearch(f,z0,pz,Hobc,EPx)       
        c1 = 1e-4;
        c2 = 0.9;
        alpha_old = 0;
        alpha_new = 1;
        alpham = 10;
        [fz0,gz0] = feval(f,z0,Hobc,EPx);
        fz_old = fz0;
        gz_old = gz0;
        i=1;
        while i<80
            z_new = z0 + alpha_new*pz;
            [fz_new,gz_new] = feval(f,z_new,Hobc,EPx);
            if (fz_new > fz0 + 2*c1*alpha_new*real(pz'*gz0)) || ((i > 1) && (fz_new >= fz_old))
                alphas = zoom(f,z0,pz,alpha_old,alpha_new,Hobc,EPx);
                return;
            end
            if abs(real(pz'*gz_new)) <= c2*abs(real(pz'*gz_old))
                alphas = alpha_new;
                return;
            end
            if real(pz'*gz_new) >= 0
                alphas = zoom(f,z0,pz,alpha_new,alpha_old,Hobc,EPx);
                return;
            end
            temp = alpha_new;
            alpha_new = 3*alpha_old + 2;
            if alpha_new > alpham
                alphas = alpham;
                return;
            end
            fz_old = fz_new;
            gz_old = gz_new;
            alpha_old = temp;
            i = i+1;
        end
        alphas = alpha_old;
    end

    function alphas = zoom(f,z0,pz,alphal,alphah,Hobc,EPx)        
        c1 = 1e-4;
        c2 = 0.9;
        [fz0,gz0] = feval(f,z0,Hobc,EPx);
        i=1;
        while i<80
            alphaj = 0.5*(alphal+alphah);
            zj = z0 + alphaj*pz;
            [fzj,gzj] = feval(f,zj,Hobc,EPx);
            fzl = feval(f,z0+alphal*pz,Hobc,EPx);
            if ((fzj > fz0 + 2*c1*alphaj*real(pz'*gz0)) || (fzj >= fzl))
                alphah = alphaj;
            else
                if abs(real(pz'*gzj)) <= c2*abs(real(pz'*gz0))
                    alphas = alphaj;
                    return;
                end
                if real(pz'*gzj)*(alphah-alphal) >= 0
                    alphah = alphal;
                end
                alphal = alphaj;
            end
            if abs(alphal-alphah)<1e-4
                alphas = alphal;
                return;
            end
            i = i+1;
        end
        alphas = 0.5*(alphal+alphah);
    end
end