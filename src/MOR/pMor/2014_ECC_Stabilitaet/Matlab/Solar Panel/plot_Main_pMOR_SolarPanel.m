%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pMOR by Matrix Interpolation
% Solar Panel
% Autor: Matthias Geuß

clear all 
clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modellerstellung und -reduktion
q = 40;
[p_ref, sys_ref_r_org, V_ref, W_ref] = Modell_SolarPanel(q);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter
p_vergl = 2;
p_int = p_ref(p_vergl);

%Auswahl = [1,3,4,5]; % liefert sehr gute Ergebnisse für Panzer + spline
Auswahl = [1,3];
p_ref = p_ref(Auswahl); % p_ref_org = [0.25 0.375 0.5 0.625 0.75]
V_ref = V_ref(Auswahl);
W_ref = W_ref(Auswahl);
sys_ref_r = sys_ref_r_org(Auswahl);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolationsmethoden

[sys_int_ri_P] = InterpolationPanzer(p_ref, p_int, V_ref, W_ref, sys_ref_r);
[sys_int_ri_S] = InterpolationPanzer_Stabil_Genau(p_ref, p_int, V_ref, W_ref, sys_ref_r);
%[sys_int_ri_SS] = InterpolationPanzer_Stabil_Genau_Spline(p_ref, p_int, V_ref, W_ref, sys_ref_r);
%[sys_int_ri_A] = InterpolationAmsallem(p_ref, p_int, V_ref, W_ref, sys_ref_r);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot
sys_vergl = sys_ref_r_org{p_vergl};

figure(1)
%omega = 10.^(2.39:0.000001:2.41);
%omega = 10.^(2:0.00001:5);
omega = 10.^(2:0.00001:5);
%omega = 10.^(2:0.000001:2.9);

sigma(sys_vergl, omega,'k','LineWidth',2);
sigma(sys_int_ri_P, omega,'k--','LineWidth',2);
%sigma(sys_int_ri_A, omega,'k-.','LineWidth',2);
%sigma(sys_int_ri_S, omega,'Color',[0.5 0.5 0.5],'LineWidth',2.5);
sigma(sys_int_ri_S, omega,'Color','r','LineWidth',2.5);
%sigma(sys_int_ri_SS, omega,'--','Color',[0.5 0.5 0.5],'LineWidth',2.5);
%xlim([10^2 2*10^4])
% xlim([10^2.9 10^2.93])
% ylim([-60 10])

size = 14; 
%set(gca, 'FontSize', size, 'FontName','Times','xtick',[10^2.39 10^2.395 10^2.4 10^2.405 10^2.41], 'xticklabel',['10^2.39', '10^2.395', '10^2.4', '10^2.405', '10^2.41']) 
set(gca, 'FontSize', size, 'FontName','Times')
%set(gca,'xtick',[10^2.39 10^2.395 10^2.4 10^2.405 10^2.41],'xticklabel',[245.5, 248.3, 251.2, 254.1, 257.0]) 
%set(gca,'xtick',[10^2.39 10^2.395 10^2.4 10^2.405 10^2.41],'xticklabel',[39.07, 39.52, 39.98, 40.44, 40.91])
h = legend('Direct ROM','Interpol. ROM ([7])','Interpol. ROM ([10])','Interpol. ROM (Geuss et al.)','Location','Southwest');
set(h, 'FontSize', size)
xlabel('Frequency (Hz)', 'FontSize', size, 'FontName','Times')
ylabel('Magnitude (dB)', 'FontSize', size, 'FontName','Times')
%xlim([10^2.39 10^2.41])
%ylim([-70 10])
grid on
box on
print -dpdf SolarPanel

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Norm

% Fehler Norm 2
e_P = (norm(sys_int_ri_P-sys_vergl,2))/norm(sys_vergl,2)
%e_A = (norm(sys_int_ri_A-sys_vergl,2))/norm(sys_vergl,2)
e_S = (norm(sys_int_ri_S-sys_vergl,2))/norm(sys_vergl,2)
%e_SS = (norm(sys_int_ri_SS-sys_vergl,2))/norm(sys_vergl,2)