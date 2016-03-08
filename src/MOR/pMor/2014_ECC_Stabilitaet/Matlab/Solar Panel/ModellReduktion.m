
function [sys_red V W] = ModellReduktion(sys, q)
% Inhalt: Reduziert das Modell 
% Autor: Matthias Geuß

% Erläuterungen:
% q: Ordnung der reduzierten Systems
% V, W: Projektionsmatrizen

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Methodenauswahl
               
%Methode = 'ICOP';
Methode = 'RK';
%Methode = 'TBR';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Berechnung

if strcmp(Methode,'ICOP')
    
    [sys_red V W aopt] = RK_ICOP(sys, [0;q]);

elseif strcmp(Methode,'RK')

    [sys_red V W] = RK(sys, [0;q] , [0;q] , @(x,y)(x'*y)); % Beidseitiger Arnoldi
    %sys_red = sss(sys_red.E\sys_red.A, sys_red.E\sys_red.B, sys_red.C);
    %[sys_red V W] = RK(sys, [0;q] , [] , @(x,y)(x'*y)); % Einseitiger Arnoldi V (W=V)
    %[sys_red V W] = RK(sys, [] , [0;q] , @(x,y)(x'*y)); % Einseitiger Arnoldi W (V=W)
    
elseif strcmp(Methode,'TBR')
    
    [sys_red V W] = TBR(sys, q); % N_Balken höchstens 5

    % [sys_red] = sss(balancmr(ss(sys)));
    % V=1;
    % W=1;
    
end

