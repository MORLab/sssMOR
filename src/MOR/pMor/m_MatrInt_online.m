function [A_int, Ap_int, B_int, C_int, E_int, x0_int, v_int, w_int] = m_MatrInt_online(R, p)
% Online Phase for matrix interpolation : Anemometer Benchmark
% by Heiko Panzer, Matthias Geuss, 2014
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ar = R{1};
Arp = R{2};
Br = R{3};
Cr = R{4};
Er = R{5};
X0 = R{6};
Vp = R{7};
Wp = R{8};
Additional_Info = R{9};
R_V = Additional_Info{1};
R_W = Additional_Info{2};
param = Additional_Info{3};
% P_lyap = Additional_Info{4};
L = Additional_Info{4};

%% Weighting functions
w = zeros(length(param),1);
w = m_Linear_weights(param, p);

[n, q] = size(Vp{1});
si = find(w); si = si(:)';

%% Reference subspace
% build V_all
tempV = zeros(n, q*length(si));
for i=1:length(si)
    tempV(:,(i-1)*q+1:i*q) = w(si(i))*Vp{si(i)}; 
end
[u,s,~] = svd(full(tempV), 'econ');
clear s v
clear R_V
R_V = full(u(:,1:q)); % first columns form R

R_W = R_V;

%% Interpolation
A_int = zeros(q,q);
B_int = zeros(q,1);
C_int = zeros(1,q);
E_int = zeros(q,q);

for i=si(:)'
  % Transformation matrix T
    [U,~,Z] = svd(Vp{i}'*R_V); 
    T = U*Z';
    
  % With M-Matrix
%     M = T'*Er{i}'*P_lyap{i};
    M = T'*Er{i}'*L{i}'*L{i};
    
    Ar_i = M*Ar{i}*T;
    Br_i = M*Br{i};
    Cr_i = Cr{i}*T;
    Er_i = M*Er{i}*T;
       
  % Interpolate
    A_int = A_int + Ar_i*w(i);
    B_int = B_int + Br_i*w(i);
    C_int = C_int + Cr_i*w(i);
    E_int = E_int + Er_i*w(i);
end

Ap_int = zeros(q,q);
x0_int = (R_V'*R_W)\R_W'*X0;
v_int = R_V'; % Definition aus taskmaster Zeile 81
w_int = R_W;

disp('H')

end
