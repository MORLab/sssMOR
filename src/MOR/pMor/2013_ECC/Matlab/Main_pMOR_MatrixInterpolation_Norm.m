%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pMOR by Matrix Interpolation
% Autor: Matthias Geuß

clear all 
clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter

p_ref = 0.4:0.4:2.4;    % Parameter der Stützstellen
p_int = 0.4:0.01:2.4;   % Parameter am Interpolationspunkt
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

[sys_int_ri_Palt] = InterpolationPanzerTest(p_ref, p_int(o), V_ref, W_ref, sys_ref_r);

[sys_int_ri_A] = InterpolationAmsallem(p_ref, p_int(o), V_ref, W_ref, sys_ref_r);

% Originalmodell / Reduziertes System
sys_org_vergl = fem_beam(p_int(o), 30);
[sys_ref_r_vergl] = ModellReduktion(sys_org_vergl, q);

% Fehler Norm 2
e_P(o) = (norm(sys_org_vergl-sys_int_ri_P,2))/norm(sys_org_vergl,2);
e_Palt(o) = (norm(sys_org_vergl-sys_int_ri_Palt,2))/norm(sys_org_vergl,2);
e_A(o) = (norm(sys_org_vergl-sys_int_ri_A,2))/norm(sys_org_vergl,2);
e_red(o) = (norm(sys_org_vergl-sys_ref_r_vergl,2))/norm(sys_org_vergl,2);
norm_org(o) = norm(sys_org_vergl,2);
norm_r_P(o) = norm(sys_int_ri_P,2);

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot
% figure(3)
% omega = 10.^(1:0.002:5);
% % bode(sys_int_o, omega,'k');
% bode(sys_int_r, omega,'k');
% % bode(sys_int_ri, omega,'r');
% % % bode(sys_ref_r_ABCE{2}, omega,'k');
% % % bode(sys_ref_r_ABCE{3}, omega,'k');
% bode(sys_int_ri_P, omega,'b');
% bode(sys_int_ri_A, omega,'r--');

figure(4)
plot(p_int,e_red,'Color',[0.5 0.5 0.5],'LineWidth',2,'MarkerSize',2,'MarkerFaceColor','k'), hold on
plot(p_int,e_Palt,'k-.','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','k'), hold on
plot(p_int,e_A,'k--','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','k'), hold on
plot(p_int,e_P,'k','LineWidth',2,'MarkerSize',2,'MarkerFaceColor','k'), hold on

grid on
%xlim([p_ref(1) p_ref(end)])
xlim([0.8 2])
ylim([0 0.25])

size = 14;

%h = legend('Panzer et al.','Amsallem et al.','Geuss et al.','Direct reduction');
h = legend('Direct ROMs','Interpolated ROMs (Panzer et al.)','Interpolated ROMs (Amsallem et al.)','Interpolated ROMs (Geuss et al.)');

set(h, 'FontSize', size, 'FontName','Times')
xlabel('Length (m)', 'FontSize', size, 'FontName','Times')
ylabel('Relative error', 'FontSize', size, 'FontName','Times')
set(gca, 'FontName','Times', 'FontSize', 14,'xtick',0.8:0.4:2.4,'ytick',0:0.05:0.3)

print -dpdf figure

