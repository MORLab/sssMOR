function [sys_int] = InterpolationManifold(p_ref, p_int, sys_ref, num_base, Methode)

% Inhalt: Interpolation im Tangentialraum 
% Autor: Matthias Geuﬂ

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter

manifold_A = 'RegulaereMatrix';
%manifold_A = 'OrthogonaleMatrix';
%manifold_A = 'SymmetrischeMatrix';

manifold_E = 'RegulaereMatrix';
%manifold_E = 'OrthogonaleMatrix';
%manifold_E = 'SymmetrischeMatrix';

manifold_B = 'ReelleMatrix';

manifold_C = 'ReelleMatrix';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Logarithmusabbildung auf den Tangentialraum

E_ref = full(sys_ref{num_base}.E);
A_ref = full(sys_ref{num_base}.A);
B_ref = full(sys_ref{num_base}.B);
C_ref = full(sys_ref{num_base}.C);

T_E_ref = {};
T_A_ref = {};
T_B_ref = {};
T_C_ref = {};
for i = 1:length(p_ref)
    
    % Abbildung der E-Matrix
    if strcmp(manifold_E,'RegulaereMatrix')
    T_E = logm(full(sys_ref{i}.E)/E_ref);
    elseif strcmp(manifold_E,'OrthogonaleMatrix')
    T_E = logm(E_ref'*full(sys_ref{i}.E));
    elseif strcmp(manifold_E,'SymmetrischeMatrix')
    T_E = logm(sqrtm(E_ref)\full(sys_ref{i}.E)/sqrtm(E_ref));
    end
    T_E = sparse(T_E);
    
    % Abbildung der A-Matrix
    if strcmp(manifold_A,'RegulaereMatrix')
    T_A = logm(full(sys_ref{i}.A)/A_ref);
    elseif strcmp(manifold_A,'OrthogonaleMatrix')
    T_A = logm(A_ref'*full(sys_ref{i}.A));
    elseif strcmp(manifold_A,'SymmetrischeMatrix')
    T_A = logm(sqrtm(A_ref)\full(sys_ref{i}.A)/sqrtm(A_ref));
    end
    T_A = sparse(T_A);
    
    % Abbildung der B-Matrix
    if strcmp(manifold_B,'ReelleMatrix')
    T_B = full(sys_ref{i}.B) - B_ref;
    end
    T_B = sparse(T_B);
    
    % Abbildung der C-Matrix
    if strcmp(manifold_C,'ReelleMatrix')
    T_C = sys_ref{i}.C - C_ref;
    end
    T_C = sparse(T_C);

%     T_E_ref = [T_E_ref {T_E}]; 
%     T_A_ref = [T_A_ref {T_A}];
%     T_B_ref = [T_B_ref {T_B}];
%     T_C_ref = [T_C_ref {T_C}];
    
    sys_ref_T{i} = sss(T_A, T_B, T_C, 0, T_E);

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolation

%Methode = 'linear';
%Methode = 'spline';
%Methode = 'pchip';

[sys_int_T] = Interpolation(p_ref, p_int, sys_ref_T, Methode);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exponentialabbildung auf die Mannigfaltigkeit

    % Abbildung der E-Matrix
    if strcmp(manifold_E,'RegulaereMatrix')
    E_int = expm(full(sys_int_T.E))*E_ref;
    elseif strcmp(manifold_E,'OrthogonaleMatrix')
    E_int = E_ref*expm(full(sys_int_T.E));
    elseif strcmp(manifold_E,'SymmetrischeMatrix')
    E_int = sqrtm(E_ref)*expm(full(sys_int_T.E))*sqrtm(E_ref);
    end
    E_int = sparse(real(E_int));

    % Abbildung der A-Matrix
    if strcmp(manifold_A,'RegulaereMatrix')
    A_int = expm(full(sys_int_T.A))*A_ref;
    elseif strcmp(manifold_A,'OrthogonaleMatrix')
    A_int = A_ref*expm(full(sys_int_T.A));
    elseif strcmp(manifold_A,'SymmetrischeMatrix')
    A_int = sqrtm(A_ref)*expm(full(sys_int_T.A))*sqrtm(A_ref);
    end
    A_int = sparse(real(A_int));
    
    % Abbildung der B-Matrix
    if strcmp(manifold_B,'ReelleMatrix')
    B_int = B_ref + full(sys_int_T.B);
    end
    B_int = sparse(B_int);
    
    % Abbildung der C-Matrix
    if strcmp(manifold_C,'ReelleMatrix')
    C_int = C_ref + full(sys_int_T.C);
    end
    C_int = sparse(C_int);
    
sys_int = sss(A_int, B_int, C_int, 0, E_int);
