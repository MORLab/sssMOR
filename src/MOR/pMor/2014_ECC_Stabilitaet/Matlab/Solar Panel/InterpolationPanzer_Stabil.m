
function [sys_int] = InterpolationPanzer_Stabil(p_ref, p_int, V, W, sys_ref)
% Inhalt: Parametrische Modellreduktion mit Matrixinterpolation, 
%         Methode von Heiko Panzer
% Autor: Matthias Geuﬂ

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bestimmung der Gewichtungsfaktoren
%w = gewichte(p_ref, p_int);
%[~,si,~] = find(w); % Eintr‰ge ungleich Null suchen

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bestimmung von Basis R

%Methode = 'Ungewichtet';
%Methode = 'Gewichtet';
Methode = 'EineBasis';

[n q] = size(V{1});

if strcmp(Methode,'Ungewichtet')
    
    tempV = zeros(n, q*length(p_ref));
    for i = 1:length(p_ref)
        tempV(:,(i-1)*q+1:i*q) = V{i};
        tempW(:,(i-1)*q+1:i*q) = W{i};
    end
    [u,~,~] = svd(full(tempV), 'econ');
    R = u(:,1:q); R = full(R);
    
    [u,~,~] = svd(full(tempW), 'econ');
    S = u(:,1:q); S=full(S);

elseif strcmp(Methode,'Gewichtet')
    
    tempV = zeros(n, q*length(si));
    for i=1:length(si)
    tempV(:,(i-1)*q+1:i*q) = w(si(i))*V{si(i)};
    %tempW(:,(i-1)*q+1:i*q) = w(si(i))*W{si(i)};
    end
    [u,~,~] = svd(full(tempV), 'econ');
    R = u(:,1:q); R = full(R);
    
%     [u,~,~] = svd(full(tempW), 'econ');
%     S = u(:,1:q); S=full(S);
    
elseif strcmp(Methode,'EineBasis')
    
    R = V{1};
        
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bestimmung von M,T und Interpolation

% A_int = zeros(q,q);
% B_int = zeros(q,1);
% C_int = zeros(1,q);
% D_int = 0;
% E_int = zeros(q,q);

for i=1:length(p_ref)
    
    [U, S, Z] = svd(V{i}'*V{1});
    Q = U*Z';

    T = (R'*V{i});
    % Minv = W{i}'*R;

    %% Trafo
    
    Ar = Q'*sys_ref{i}.A*Q;
    Br = Q'*sys_ref{i}.B;
    Cr = sys_ref{i}.C*Q;
    
    
%     Ar = T*sys_ref{i}.A/T;
%     Br = T*sys_ref{i}.B;
%     Cr = sys_ref{i}.C/T;

    %% Stabilisierung
    
    RS = eye(length(Ar));
    P = lyap(Ar',RS); 
    %L = chol(P);
     L = lyapchol(Ar',RS);
    %L = sqrtm(P);   
    
    Ar_Stab = L*Ar/L;
    Br_Stab = L*Br;
    Cr_Stab = Cr/L;
    sys_ref_Kon{i} = sss(Ar_Stab, Br_Stab, Cr_Stab);

%     Ar_Stab = P*Ar;
%     Br_Stab = P*Br;
%     Cr_Stab = Cr;
%     Er_Stab = P;
%     sys_ref_Kon{i} = sss(Ar_Stab, Br_Stab, Cr_Stab, 0 , Er_Stab);

%     Ar_Stab = P*Ar/P;
%     Br_Stab = P*Br;
%     Cr_Stab = Cr/P;
%     sys_ref_Kon{i} = sss(Ar_Stab, Br_Stab, Cr_Stab);

%     Ar_Stab = T*L*sys_ref{i}.A/L/T;
%     Br_Stab = T*L*sys_ref{i}.B;
%     Cr_Stab = sys_ref{i}.C/L/T;
%     sys_ref_Kon{i} = sss(Ar_Stab, Br_Stab, Cr_Stab);
    
end

[sys_int] = Interpolation(p_ref, p_int, sys_ref_Kon, 'linear');
%[sys_int] = Interpolation(p_ref, p_int, sys_ref_Kon, 'spline');
%[sys_int] = Interpolation(p_ref, p_int, sys_ref_Kon, 'pchip');
%[sys_int] = InterpolationManifold(p_ref, p_int, sys_ref_Kon, 1, 'spline');


