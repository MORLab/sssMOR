%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pMOR by Matrix Interpolation
% Autor: Matthias Geuß

clear all 
clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter
% p_int_vec = 0.8:0.01:1.2;
% p_ref = 0.8:0.2:1.2;    % Parameter der Stützstellen
%p_int_vec = 0.8:0.002:1.2;
p_int_vec = 0.8:0.01:1.2;
p_ref = 0.8:0.1:1.2;    % Parameter der Stützstellen
q = 23;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modellerstellung der Originalmodelle

% Originalsysteme an Stützstellen
sys_ref_o = {};
for i = 1:length(p_ref)
    
    sys = ModellErstellung(p_ref(i));
    sys_ref_o = [sys_ref_o {sys}];
    clear sys;
  
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modellreduktion der lokalen Modelle

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

for ip = 1: length(p_int_vec)
p_int = p_int_vec(ip);              % Parameter am Interpolationspunkt

% Originalsystem am Interpolationspunkt
[sys_int_o] = ModellErstellung(p_int);
% Reduziertes Modell am Interpolationspunkt
[sys_int_r] = ModellReduktion(sys_int_o, q);
   
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolationsmethoden
%[sys_int_ri_P] = InterpolationPanzer_PvorQ(p_ref, p_int, V_ref, W_ref, sys_ref_r);
[sys_int_ri_P] = InterpolationPanzer(p_ref, p_int, V_ref, W_ref, sys_ref_r);
[sys_int_ri_S] = InterpolationPanzer_Stabil_Genau(p_ref, p_int, V_ref, W_ref, sys_ref_r);
[sys_int_ri_A] = InterpolationAmsallem(p_ref, p_int, V_ref, W_ref, sys_ref_r);

%plot_eig_A(ip) = max(real(eig(sys_int_ri_A)));
plot_eig_P(ip) = max(real(eig(sys_int_ri_P)));
plot_eig_S(ip) = max(real(eig(sys_int_ri_S)));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot

% figure(1)
% omega = 10.^(0:0.002:6);
% sigma(sys_int_r, omega,'k');
% sigma(sys_int_ri_P, omega,'r--');
% sigma(sys_int_ri_A, omega,'b--');
% sigma(sys_int_ri_S, omega,'g--');
% ylim([-100 40])

sys_org_vergl = fem_beam(p_int_vec(ip), 12);
[sys_ref_r_vergl] = ModellReduktion(sys_org_vergl, q);
%e_S(ip) = (norm(sys_org_vergl-sys_int_ri_S,2))/norm(sys_org_vergl,2);
e_S(ip) = (norm(sys_ref_r_vergl-sys_int_ri_S,2))/norm(sys_ref_r_vergl,2);
e_P(ip) = (norm(sys_ref_r_vergl-sys_int_ri_P,2))/norm(sys_ref_r_vergl,2);
e_A(ip) = (norm(sys_ref_r_vergl-sys_int_ri_A,2))/norm(sys_ref_r_vergl,2);


end

figure(8)
plot(p_int_vec,plot_eig_S,'k','LineWidth',2); hold on
plot(p_int_vec,plot_eig_P,'k--','LineWidth',2); hold on
%plot(p_int_vec,plot_eig_A,'k--','LineWidth',2); hold on

%plot([0.8 1.2],[0 0],'k--');

size = 14; 
%set(gca, 'FontSize', size, 'FontName','Times','xtick',[10^2.39 10^2.395 10^2.4 10^2.405 10^2.41], 'xticklabel',['10^2.39', '10^2.395', '10^2.4', '10^2.405', '10^2.41']) 
set(gca, 'FontSize', size, 'FontName','Times')
%set(gca,'xtick',[10^2.39 10^2.395 10^2.4 10^2.405 10^2.41],'xticklabel',[245.5, 248.3, 251.2, 254.1, 257.0]) 
h = legend('Interpolated ROMs (Geuss et al.)','Interpolated ROMs (Reference)','Location','NorthEast');
set(h, 'FontSize', size)
xlabel('Length (m)', 'FontSize', size, 'FontName','Times')
ylabel('Real part of the rightmost eigenvalue', 'FontSize', size, 'FontName','Times')
xlim([0.9 1.1])
ylim([-4.5 0.5])
grid on
box on
print -dpdf Beam2

figure(9)
hold on; plot(p_int_vec,e_S,'k','LineWidth',1.5,'MarkerSize',2,'MarkerFaceColor','k'), hold on
plot(p_int_vec,e_P,'k--','LineWidth',1.5,'MarkerSize',2,'MarkerFaceColor','k'), hold on
%plot(p_int_vec,e_A,'k--','LineWidth',1.5,'MarkerSize',2,'MarkerFaceColor','g')
%plot(p_int_vec,e_red,'k','LineWidth',1.5,'MarkerSize',2,'MarkerFaceColor','k')
size = 14; 
%set(gca, 'FontSize', size, 'FontName','Times','xtick',[10^2.39 10^2.395 10^2.4 10^2.405 10^2.41], 'xticklabel',['10^2.39', '10^2.395', '10^2.4', '10^2.405', '10^2.41']) 
set(gca, 'FontSize', size, 'FontName','Times')
%set(gca,'xtick',[10^2.39 10^2.395 10^2.4 10^2.405 10^2.41],'xticklabel',[245.5, 248.3, 251.2, 254.1, 257.0]) 
h = legend('Interpolated ROMs (Geuss et al.)', 'Interpolated ROMs (Reference)','Location','NorthEast');
% set(h, 'FontSize', size)
xlabel('Length (m)', 'FontSize', size, 'FontName','Times')
ylabel('Relative error', 'FontSize', size, 'FontName','Times')
xlim([0.9 1.1])
ylim([0 3])
grid on
box on
print -dpdf Beam1