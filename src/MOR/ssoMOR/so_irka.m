function [sysr, V, iter] = so_irka(sys,nr,Opts)
% Iterative-Rational-Krylov-Algorithmus f¨¹r
% Systeme zweiter Ordnung (SO-IRKA)
%
% Input    -      tol:     Tolernz der Hauptschleife   
%                  nr:     Reduktionsordnung
%                nmax:     Maximale Anzahl der Iterationen
%                   M:     Massenmatrix des Originalsystems
%                   K:     Steifigkeitsmatrix des Originalsystems
%                   D:     D?mpfungsmatrix des Originalsystems
%                   B:     Eingangsmatrix des Originalsystems
%
% 
% Output    -       V:     Rechte Projektionsmatrix
%                sysr:     Reduziertes Modell
%                iter:     Anzahl der Iterationen

%%  Define optional execution parameters
Def.tol     = 1e-3;
Def.nmax    = 50;

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end  

%% Initialize
[M,D,K,B]   = dssdata(sys);
V_ini       = [speye(nr);sparse(sys.n-nr,nr)];
sysr        = projectiveMor(sys,V_ini);
[M_ini,D_ini,K_ini,B_ini] = dssdata(sysr);

%% Run SO-IRKA iteration
numfree = size(B,1);
n = 0;
[EP,Lb_n] = Neueptv(full(M_ini),full(D_ini),full(K_ini),full(B_ini));
V = Neuvr(EP,Lb_n);
es = 1;
while n <=Opts.nmax && es > Opts.tol
    M_n = V'*M*V;
    K_n = V'*K*V;
    D_n = V'*D*V;
    B_n = V'*B; 
    [EP_n,Lb_n] = Neueptv(M_n,D_n,K_n,B_n);
    V = Neuvr(EP_n,Lb_n);
    es = norm(EP_n-EP)/norm(EP); 
    EP = EP_n;
    n = n+1;
end

iter = n;
sysr = projectiveMor(sys,V);

%% ================================================================
%       AUXILIARY
%  ================================================================

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
        V_n = zeros(numfree,nr);
        i = 1;
        while i <= nr
            K_sigma = EP_n(i)^2*M+EP_n(i)*D+K;
            [L,U,P,Q,R] = lu(K_sigma);
            V_temp = Q*(U\(L\(P*(R\(B*Lb_n(:,i))))));
            if ~isreal(EP_n(i))
                V_n(:,i:i+1) = [real(V_temp) imag(V_temp)];
                i = i + 2;
            else
                V_n(:,i) = V_temp;
                i = i + 1;
            end
%             [V_rn(:,1:i-1),~] = qr(V_rn(:,1:i-1),0);
        end
        [V_n,~] = qr(V_n,0);
    end
end