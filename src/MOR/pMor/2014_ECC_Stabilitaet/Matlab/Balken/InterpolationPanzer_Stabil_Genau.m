
function [sys_int] = InterpolationPanzer_Stabil_Genau(p_ref, p_int, V, W, sys_ref)
% Inhalt: Parametrische Modellreduktion mit Matrixinterpolation, 
%         Methode von Heiko Panzer
% Autor: Matthias Geuß

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bestimmung von M,T und Interpolation
L={};
for i=1:length(p_ref)
    
    Ar = sys_ref{i}.A;
    Br = sys_ref{i}.B;
    Cr = sys_ref{i}.C;     
    
    %% Stabilisierung
   %RS = Cr'*Cr;
   %RS = (V{i}'*W{2})*(V{i}'*W{2})';
   %RS = V{i}'*W{2}*W{2}'*V{i};
    RS = eye(length(Ar));
      %    Ra = rand(length(Ar));
    %    RS = Ra'*Ra;
    P = lyap(Ar',RS);
    
% A = full(Ar);
% Ps = sdpvar(length(A),length(A));
% F = [Ps >= 0, A'*Ps+Ps*A <= -eye(length(A),length(A))*0.001];
% %F = [Ps >= 0, Ps*A+A'*Ps <= 0, ];
% %solvesdp(F,norm(V{2}*Ps-V_neu{1},'fro'));
% solvesdp(F,trace(Ps))
% trace(Ps)
% %solvesdp(F,norm(V{2}-V_neu{1}*Ps,'fro'));
% %solvesdp(F,if double(Ps)== NaN, norm(V{2}-V_neu{1},'fro'); else norm(V{2}*sqrtm(double(Ps))-V_neu{1},'fro');end);
% P = double(Ps);
    
    
%     
% A = full(Ar);
% k = length(A);
% cvx_begin sdp
% cvx_solver sdpt3
% cvx_solver SeDuMi
% variable P(k,k) symmetric
% minimize(trace(P));
% subject to
% P - 10^-8*eye(k,k)>= 0;
% P*A+A'*P <= 0;
% cvx_end
    
%    L{i} = chol(P);
    %L{i} = eye(23,23);
    L{i} = lyapchol(Ar',RS);
    %L{i} = sqrtm(P);
    
%     Ar_Stab = L*Ar/L;
%     Br_Stab = L*Br;
%     Cr_Stab = Cr/L;
%     sys_ref_Stab{i} = sss(Ar_Stab, Br_Stab, Cr_Stab);
    
    V_neu{i} = V{i}/L{i};

end


for i=1:length(p_ref)
    
    [U, S, Z] = svd(V_neu{i}'*V_neu{2});
    Q = U*Z';
    J_alt = norm(V_neu{i} - V_neu{2},'fro')
    J = norm(V_neu{i}*Q - V_neu{2},'fro')

%Q = eye(length(Ar),length(Ar));

    %% Trafo
    
    Ar_Genau = Q'*L{i}*sys_ref{i}.A/L{i}*Q;
    Br_Genau = Q'*L{i}*sys_ref{i}.B;
    Cr_Genau = sys_ref{i}.C/L{i}*Q; 
    
%     Ar_Genau = Q'*sys_ref_Stab{i}.A*Q;
%     Br_Genau = Q'*sys_ref_Stab{i}.B;
%     Cr_Genau = sys_ref_Stab{i}.C*Q;  
    
    sys_ref_Kon{i} = sss(Ar_Genau, Br_Genau, Cr_Genau);
    
end

[sys_int] = Interpolation(p_ref, p_int, sys_ref_Kon, 'linear');
%[sys_int] = Interpolation(p_ref, p_int, sys_ref_Kon, 'spline');
%[sys_int] = Interpolation(p_ref, p_int, sys_ref_Kon, 'pchip');
%[sys_int] = InterpolationManifold(p_ref, p_int, sys_ref_Kon, 1, 'spline');


