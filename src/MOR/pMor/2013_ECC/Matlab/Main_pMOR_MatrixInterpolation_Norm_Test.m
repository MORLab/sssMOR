%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pMOR by Matrix Interpolation
% Autor: Matthias Geuß

clear all 
clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter

p_ref = 0.4:0.4:2.4;    % Parameter der Stützstellen
p_int = 0.4:0.05:2.4;   % Parameter am Interpolationspunkt
q = 10;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modellerstellung der Originalmodelle

% Originalsystem am Interpolationspunkt
%[sys_int_o] = ModellErstellung(p_int);

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
%[sys_int_r] = ModellReduktion(sys_int_o);

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
    
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolationsmethoden

for o = 1:1:length(p_int)

[sys_int_ri_P] = InterpolationPanzer(p_ref, p_int(o), V_ref, W_ref, sys_ref_r);

%[sys_int_ri_P] = InterpolationPanzerTest(p_ref, p_int(o), V_ref, W_ref, sys_ref_r);

[sys_int_ri_A] = InterpolationAmsallem(p_ref, p_int(o), V_ref, W_ref, sys_ref_r);

% Originalmodell / Reduziertes System
sys_org_vergl = fem_beam(p_int(o), 10);
[sys_ref_r_vergl] = ModellReduktion(sys_org_vergl, q);

% Fehler Norm 2
e_P(o) = (norm(sys_int_ri_P-sys_org_vergl,2))/norm(sys_org_vergl,2);
e_A(o) = (norm(sys_int_ri_A-sys_org_vergl,2))/norm(sys_org_vergl,2);
e_red(o) = (norm(sys_ref_r_vergl-sys_org_vergl,2))/norm(sys_org_vergl,2);
norm_org(o) = norm(sys_org_vergl,2);
norm_r_P(o) = norm(sys_int_ri_P,2);

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot

figure(1)
plot(p_int,e_P,'b-o','MarkerSize',2,'MarkerFaceColor','b'), hold on
plot(p_int,e_A,'r--o','MarkerSize',2,'MarkerFaceColor','r')
plot(p_int,e_red,'g-o','MarkerSize',2,'MarkerFaceColor','g')
grid on
xlim([p_ref(1) p_ref(end)])
ylim([0 1.2])

