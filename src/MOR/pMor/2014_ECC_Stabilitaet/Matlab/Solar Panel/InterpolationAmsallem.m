function [sys_int] = InterpolationAmsallem(p_ref, p_int, V, W, sys_ref)
% Inhalt: Parametrische Modellreduktion mit Matrixinterpolation, 
%         Methode von David Amsallem
% Autor: Matthias Geuﬂ


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Umrechnung der Systeme in ABC statt ABCE

for i = 1:length(p_ref)
    sys_ref_ABC{i} = sss(sys_ref{i}.E\sys_ref{i}.A, sys_ref{i}.E\sys_ref{i}.B, sys_ref{i}.C);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basissystem aus p_ref

num_base = 1;
V_base = V{num_base};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kongruenztransformation

Q = {};
for i = 1:length(p_ref)

    [U, S, Z] = svd(V{i}'*V_base);
    Q{i} = U*Z';
    sys_ref_A = Q{i}'*sys_ref_ABC{i}.A*Q{i};
    sys_ref_B = Q{i}'*sys_ref_ABC{i}.B;
    sys_ref_C = sys_ref_ABC{i}.C*Q{i};
    
    sys_ref_Kon{i} = sss(sys_ref_A, sys_ref_B, sys_ref_C);

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolation auf Mannigfaltigkeit

%[sys_int] = InterpolationManifold(p_ref, p_int, sys_ref_Kon, num_base, 'linear');
[sys_int] = InterpolationManifold(p_ref, p_int, sys_ref_Kon, num_base, 'spline');