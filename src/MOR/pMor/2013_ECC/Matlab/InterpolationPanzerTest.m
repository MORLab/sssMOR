
function [sys_int] = InterpolationPanzerTest(p_ref, p_int, V, W, sys_ref)
% Inhalt: Parametrische Modellreduktion mit Matrixinterpolation, 
%         Methode von Heiko Panzer
% Autor: Matthias Geuﬂ

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bestimmung der Gewichtungsfaktoren
%w = gewichte(p_ref, p_int);
%[~,si,~] = find(w); % Eintr‰ge ungleich Null suchen

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bestimmung von Basis R

Methode = 'Ungewichtet';
%Methode = 'Gewichtet';
%Methode = 'EineBasis';

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
    
%      [u,~,~] = svd(full(tempW), 'econ');
%      S = u(:,1:q); S=full(S);
    
elseif strcmp(Methode,'EineBasis')
    
    R = V{3};
        
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bestimmung von M,T und Interpolation

% A_int = zeros(q,q);
% B_int = zeros(q,1);
% C_int = zeros(1,q);
% D_int = 0;
% E_int = zeros(q,q);

for i=1:length(p_ref)
    
%     [U, S, Z] = svd(V{i}'*V{1});
%     Q = U*Z';
%     [U2, S2, Z2] = svd(W{i}'*W{1});
%     Q2 = U2*Z2';

    T = (R'*V{i});
    %Minv = W{i}'*R;
    Minv = W{i}'*R;
    
   % Minv = W{i}'*W{1};   
    
    Ar = Minv\sys_ref{i}.A/T;
    Br = Minv\sys_ref{i}.B;
    Cr = sys_ref{i}.C/T;
    Er = Minv\sys_ref{i}.E/T;
    
%     Ar = Minv\sys_ref{i}.A*Q;
%     Br = Minv\sys_ref{i}.B;
%     Cr = sys_ref{i}.C*Q;
%     Er = Minv\sys_ref{i}.E*Q;

%     Ar = sys_ref{i}.A/T;
%     Br = sys_ref{i}.B;
%     Cr = sys_ref{i}.C/T;
%     Er = sys_ref{i}.E/T;

%     Ar = sys_ref{i}.A*Q;
%     Br = sys_ref{i}.B;
%     Cr = sys_ref{i}.C*Q;
%     Er = sys_ref{i}.E*Q;

%     EinvA = sys_ref{i}.E\sys_ref{i}.A;
%     EinvB = sys_ref{i}.E\sys_ref{i}.B;
%     Ar = Q'*EinvA*Q;
%     Br = Q'*EinvB;
%     Cr = sys_ref{i}.C*Q;
%     Er = eye(q);

%     EinvA = sys_ref{i}.E\sys_ref{i}.A;
%     EinvB = sys_ref{i}.E\sys_ref{i}.B;
%     Ar = T*EinvA/T;
%     Br = T*EinvB;
%     Cr = sys_ref{i}.C/T;
%     Er = eye(q);

% 
%     EinvA = sys_ref{i}.E\sys_ref{i}.A;
%     EinvB = sys_ref{i}.E\sys_ref{i}.B;
%     Ar = T\EinvA*T;
%     Br = T\EinvB;
%     Cr = sys_ref{i}.C*T;
%     Er = eye(q);

%     Ar = Q2*sys_ref{i}.A*Q;
%     Br = Q2*sys_ref{i}.B;
%     Cr = sys_ref{i}.C*Q;
%     Er = Q2*sys_ref{i}.E*Q;
    
%     Ar = Minv\sys_ref{i}.A*Q;
%     Br = Minv\sys_ref{i}.B;
%     Cr = sys_ref{i}.C*Q;
%     Er = Minv\sys_ref{i}.E*Q;

    
%     % Interpolation
%     A_int = A_int + Ar*w(i);
%     B_int = B_int + Br*w(i);
%     C_int = C_int + Cr*w(i);
%     E_int = E_int + Er*w(i); 
    
    sys_ref_Kon{i} = sss(Ar, Br, Cr, 0, Er);
    
end
%sys_int = sss(A_int, B_int, C_int, 0 , E_int);

[sys_int] = Interpolation(p_ref, p_int, sys_ref_Kon, 'linear');
%[sys_int] = Interpolation(p_ref, p_int, sys_ref_Kon, 'spline');
%[sys_int] = Interpolation(p_ref, p_int, sys_ref_Kon, 'pchip');
%[sys_int] = InterpolationManifold(p_ref, p_int, sys_ref_Kon, 1, 'spline');


