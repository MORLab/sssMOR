function R = m_MatrInt_offline(A,Ap,B,C,D,E,X0,U,P,h,T,param,r)
% Offline Phase for matrix interpolation: 
% by Heiko Panzer, Matthias Geuss, 2014
% Aenderungen, Erweiterungen und Korrekturen (Maria Cruz Varona 2015):
%- Two-Step-Approach, One-Sided Rational Krylov hinzugefügt
%- Zeilen-Dimension von R_V und R_W angepasst 

% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sample Points
anz_s = length(param);

% Definitons
Vp = cell(1,length(param));
Wp = cell(1,length(param));
Ar = cell(1,length(param));
Br = cell(1,length(param));
Cr = cell(1,length(param));
Dr = cell(1,length(param));
Er = cell(1,length(param));
Apr = 0;
% P_lyap = cell(1,length(param));
L = cell(1,length(param));

n = length(A);
Vall = zeros(n, r*anz_s);

for j=1:anz_s
       
    sys = sss((A+param(j)*Ap),B,C,0,E); 
    
%% Local Reduction
    %% Twosided IRKA
      % initial values for IRKA
    s0 = -eigs(A+param(j)*Ap,E,r,0, struct('tol',1e-6,'v0',sum(B,2)));
    try
        s0 = cplxpair(s0);
    catch e
        s0(imag(s0)~=0) = 0;
    end
    
    % initial values for IRKA
    s0 = -eigs(A+param(j)*Ap,E,r,0, struct('tol',1e-6,'v0',sum(B,2)));
    if mod(nnz(imag(s0)),2)~=0  % complex conjugate
        % find the s0 which has no compl.conj. partner
        s0(imag(s0) == imag(sum(s0))) = 0;
    end
    s0= cplxpair(s0);


    [~,V,W] = IRKA(sys, s0, 50, 1e-7);
    Vp{j} = V;
    Wp{j} = W;

    %% Two-Step-Approach (Krylov + BT):
    %1.) One-Sided Krylov
%     r1 = 100; % Intermediate reduced order, anemometer: 100, synthetic: 200 / 250
%     s0 = 0; % Shift
%     [sysr1, V1, W1] = RK(sys, [s0;r1], [], @(x,y)(x'*y)); % One-sided Arnoldi V (W=V)
%     check_SD(sysr1);
%     
%     %2.) Balanced Truncation
%     [sysr2, V2, W2] = TBR(sysr1, r);
% %     [sysr2,red_info] = balancmr(ss(sysr1),r); 
% %     [sysr2, g] = balreal(ss(sysr1))
% 
%     %Assembly
%     Vp{j} = orth(V1*V2);
%     Wp{j} = orth(W1*W2);
%     
%     clear V1 W1 V2 W2 sysr1 sysr2
    
    %% One-Sided Rational Krylov:
%     s0 = 0; % Shift
%     [~, V, W] = RK(sys, [s0;r], [], @(x,y)(x'*y)); % One-sided Arnoldi V (W=V)
% 
%     Vp{j} = V;
%     Wp{j} = W;
    
%% Model reduction  
    Ar{j} = Wp{j}'*(A+param(j)*Ap)*Vp{j};
    Br{j} = Wp{j}'*B;
    Cr{j} = C*Vp{j};
    Dr{j} = 0;
    Er{j} = Wp{j}'*E*Vp{j}; 
    
%% Lyapunov-Equations
    %% lyap
%     P_lyap{j} = lyap(Ar{j}', Cr{j}'*Cr{j}, [], Er{j}');
%     P_lyap{j} = lyap(Ar{j}', eye(r), [], Er{j}');

    %% lyapchol
    % Liefert eine Fehlermeldung, wenn ROM nicht stabil.
%     L{j} = lyapchol(Ar{j}', Cr{j}', Er{j}');
%     L{j} = lyapchol(Ar{j}', eye(r), Er{j}');
    L{j} = lyapchol(Ar{j}', diag(Cr{j})/norm(Cr{j}), Er{j}');
    
%% Calculation of Q_all
    Vall(:,(j-1)*r+1:j*r) = Vp{j};

%     figure
%     sigma(sss(Ar{j},Br{j},Cr{j},0,Er{j}), 10.^(-2:0.001:4));
%     hold on
%     sigma(sys, 10.^(0:0.05:4), 'k*');
%     drawnow

end

% Falls R_V, R_W offline definiert
R_V = zeros(length(A),r);
R_W = zeros(length(A),r);

Q_all = Vall'*Vall;

% Additional_Info = {R_V R_W param P_lyap Q_all};
Additional_Info = {R_V R_W param L Q_all};

R = {Ar, Apr, Br, Cr, Er, X0, Vp, Wp, Additional_Info, @m_MatrInt_online_HP2};
