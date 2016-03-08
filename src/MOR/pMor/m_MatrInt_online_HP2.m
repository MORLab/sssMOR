function [A_int, Ap_int, B_int, C_int, E_int, x0_int, v_int, w_int] = m_MatrInt_online_HP2(R, p)
% Online Phase for matrix interpolation : Anemometer Benchmark
% by Heiko Panzer, Matthias Geuss, 2014 - Maria Cruz Varona 2015
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ar = R{1};
% Arp = R{2};
Br = R{3};
Cr = R{4};
Er = R{5};
X0 = R{6};
Vp = R{7};
% Wp = R{8};
Additional_Info = R{9};
% R_V = Additional_Info{1};
% R_W = Additional_Info{2};
param = Additional_Info{3};
% P_lyap = Additional_Info{4};
L = Additional_Info{4};
Q_all = Additional_Info{5};

%% Weighting functions
% w = zeros(length(param),1);
w = m_Linear_weights(param, p)

%%
[n, q] = size(Vp{1});
si = find(w); si = si(:)';

R = zeros(n,q);
% R = [];

wall = kron(diag(w), eye(q));
wall(:,diag(wall)==0)=[];
[~,s,z] = svd(wall'*Q_all*wall); s = sqrt(s);
if length(s)==q
    dissim = 0;
else
    dissim = (s(q+1,q+1)/s(1,1))*100
end

%% Interpolation
A_int = zeros(q,q);
B_int = zeros(q,1);
C_int = zeros(1,q);
E_int = zeros(q,q);

temp = z/s; % for calculation of R
% 
% RHS = zeros(1,q);
% for i=length(si)
%     RHS = RHS + w(si(i)) * Cr{si(i)};
% end

for i=1:length(si)
  % Transformation matrix T
    RtV = s(1:q,:)*z((i-1)*q+1:i*q,:)'/w(si(i));
%     [U,~,Z] = svd(RtV');
%     T = U*Z';
    T = RtV';
    
  % With M-Matrix
%     M = T'*Er{si(i)}'*P_lyap{si(i)};
%     M = T'*Er{i}'*L{si(i)}'*L{si(i)};

%     L = lyapchol(Ar{si(i)}', RHS', Er{si(i)}');
%     Lx = lyapchol(Ar{si(i)}', diag(Cr{si(i)}), Er{si(i)}');
    Lx = L{si(i)};
    
    %     L = lyapchol(Ar{si(i)}', Cr{si(i)}', Er{si(i)}');
%     L = lyapchol(Ar{si(i)}',T, Er{si(i)}');
%     Px = lyap(Ar{si(i)}', RHS'*RHS+0.01*eye(q), Er{si(i)}');
%     M = T'*Er{si(i)}'*Px;

% %     Ar_i = M*Ar{si(i)}*T;
%     Ar_i = T'*((L{si(i)}*Er{si(i)})'*(L{si(i)}*Ar{si(i)}))*T;
% %     Br_i = M*Br{si(i)};
%     Br_i = T'*(Er{si(i)}'*(L{si(i)}'*(L{si(i)}*Br{si(i)})));
%     Cr_i = Cr{si(i)}*T;
% %     Er_i = M*Er{si(i)}*T;
%     Er_i = (L{si(i)}*Er{si(i)}*T)'*(L{si(i)}*Er{si(i)}*T);

    Ar_i = T'*((Lx*Er{si(i)})'*(Lx*Ar{si(i)}))*T;
    Br_i = T'*(Er{si(i)}'*(Lx'*(Lx*Br{si(i)})));
    Cr_i = Cr{si(i)}*T;
    Er_i = (Lx*Er{si(i)}*T)'*(Lx*Er{si(i)}*T);
    
  % Interpolate
    A_int = A_int + Ar_i*w(si(i));
    B_int = B_int + Br_i*w(si(i));
    C_int = C_int + Cr_i*w(si(i));
    E_int = E_int + Er_i*w(si(i));
    
  % Compute R
    R = R + Vp{si(i)}*temp((i-1)*q+1:i*q,1:q)*w(si(i));
end
% R_test = Vall*wall*z/s;
R_V = R;
R_W = R_V;

Ap_int = zeros(q,q);
x0_int = (R_V'*R_W)\R_W'*X0;
v_int = R_V'; % Definition aus taskmaster Zeile 81
w_int = R_W;
% 
% figure
% sigma(sss(A_int,B_int,C_int,0,E_int), 10.^(-1:0.01:4));
% hold on

end
