
function [sys_int] = InterpolationPanzer_Stabil_Genau_Spline(p_ref, p_int, V, W, sys_ref)
% Inhalt: Parametrische Modellreduktion mit Matrixinterpolation, 
%         Methode von Heiko Panzer
% Autor: Matthias Geuﬂ

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bestimmung von M,T und Interpolation

for i = 1:length(p_ref)
    sys_ref_ABC{i} = sss(sys_ref{i}.E\sys_ref{i}.A, sys_ref{i}.E\sys_ref{i}.B, sys_ref{i}.C);
end

L={};
for i=1:length(p_ref)
    
    Ar = sys_ref_ABC{i}.A;
    Br = sys_ref_ABC{i}.B;
    Cr = sys_ref_ABC{i}.C;     
    
    %% Stabilisierung
    
    RS = eye(length(Ar));
%     Ra = rand(length(Ar));
%     RS = Ra'*Ra;
    P = lyap(Ar',RS);
    L{i} = chol(P);
    %L = lyapchol(Ar',RS);
    %L{i} = sqrtm(P);
    
%     Ar_Stab = L*Ar/L;
%     Br_Stab = L*Br;
%     Cr_Stab = Cr/L;
%     sys_ref_Stab{i} = sss(Ar_Stab, Br_Stab, Cr_Stab);
    
    V_neu{i} = V{i}/L{i};

end


for i=1:length(p_ref)
    
    [U, S, Z] = svd(V_neu{i}'*V_neu{1});
    Q = U*Z';
    J_alt = norm(V_neu{i} - V_neu{1},'fro')
    J = norm(V_neu{i}*Q - V_neu{1},'fro')

    %% Trafo
    
    Ar_Genau = Q'*L{i}*sys_ref_ABC{i}.A/L{i}*Q;
    Br_Genau = Q'*L{i}*sys_ref_ABC{i}.B;
    Cr_Genau = sys_ref_ABC{i}.C/L{i}*Q; 
    
%     Ar_Genau = Q'*sys_ref_Stab{i}.A*Q;
%     Br_Genau = Q'*sys_ref_Stab{i}.B;
%     Cr_Genau = sys_ref_Stab{i}.C*Q;  
    
    sys_ref_Kon{i} = sss(Ar_Genau, Br_Genau, Cr_Genau);
    
end

%[sys_int] = Interpolation(p_ref, p_int, sys_ref_Kon, 'linear');
[sys_int] = Interpolation(p_ref, p_int, sys_ref_Kon, 'spline');
%[sys_int] = Interpolation(p_ref, p_int, sys_ref_Kon, 'pchip');
%[sys_int] = InterpolationManifold(p_ref, p_int, sys_ref_Kon, 1, 'spline');


