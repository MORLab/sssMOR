%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pMOR by Matrix Interpolation
% Autor: Matthias Geuß

clear all 
clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter

p_ref = 0.8:0.4:1.2;    % Parameter der Stützstellen
p_int = 1;              % Parameter am Interpolationspunkt
q = 10;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modellerstellung der Originalmodelle

% Originalsystem am Interpolationspunkt
[sys_int_o] = ModellErstellung(p_int);

% Originalsysteme an Stützstellen
sys_ref_o = {};
for i = 1:length(p_ref)
    
    sys = ModellErstellung(p_ref(i));
    sys_ref_o = [sys_ref_o {sys}];
    clear sys;
  
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modellreduktion der lokalen Modelle

% Reduziertes Modell am Interpolationspunkt
[sys_int_r] = ModellReduktion(sys_int_o, q);

% Reduzierte Modelle an Stützstellen
V_ref = {};
W_ref = {};
sys_ref_r = {};
for i = 1:length(p_ref)
    
    [sys_red V W] = ModellReduktion(sys_ref_o{i}, q);
    V_ref = [V_ref {V}];
    W_ref = [W_ref {W}];
    sys_ref_r = [sys_ref_r {sys_red}];
    clear V W sys_red;
    
    temp = max(real(eig(sys_ref_r{1,i})));
    
    if (temp>0)
        i %check if a reduced model is not asympt. stable
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolationsmethoden

[sys_int_ri_P] = InterpolationPanzer(p_ref, p_int, V_ref, W_ref, sys_ref_r);
[sys_int_ri_S] = InterpolationPanzer_Stabil_Genau(p_ref, p_int, V_ref, W_ref, sys_ref_r);
[sys_int_ri_A] = InterpolationAmsallem(p_ref, p_int, V_ref, W_ref, sys_ref_r);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot

figure(1)
omega = 10.^(0:0.002:6);
sigma(sys_int_r, omega,'k');
sigma(sys_int_ri_P, omega,'r');
sigma(sys_int_ri_A, omega,'y--');
sigma(sys_int_ri_S, omega,'g--');
ylim([-160 -30])

e_P = (norm(sys_int_ri_P-sys_int_r,2))/norm(sys_int_r,2)
e_A = (norm(sys_int_ri_A-sys_int_r,2))/norm(sys_int_r,2)
e_S = (norm(sys_int_ri_S-sys_int_r,2))/norm(sys_int_r,2)