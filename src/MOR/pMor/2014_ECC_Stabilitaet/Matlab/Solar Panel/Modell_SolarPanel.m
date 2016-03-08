function [p_ref, sys_ref_r_org, V_ref, W_ref] = Modell_SolarPanel(q)

% Inhalt: Reduziert das Modell "Solar Panel" für verschiedene 
%         "thickness of facesheet [mm]" von p = [0.25 0.375 0.5 0.625 0.75]
% Autor: Matthias Geuß

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition

p_ref = [0.25 0.375 0.5 0.625 0.75];
V_ref = {};
W_ref = {};
sys_ref_r_org = {};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reduktion


% p = 0.25
cd ModelsSolarPanel
load sys3_1.mat
cd ..

M_sys = sys.M;
K_sys = sys.K;
%D_sys = sys.D;
%D_sys = 10^(-8)*(M_sys+K_sys);
D_sys = 0.000001*K_sys + 0.000001*M_sys;
B_sys = sys.B;
C_sys = sys.Cp;

o = length(M_sys);

E = [K_sys, zeros(o,o); zeros(o,o), M_sys];
A = [zeros(o,o), K_sys; -K_sys, -D_sys];
B = [zeros(o,1); B_sys];
C = [C_sys , zeros(1,o)];

sys_org = sss(A,B,C,0,E);
[sys_red V W] = ModellReduktion(sys_org,q);
V_ref = [V_ref {V}];
W_ref = [W_ref {W}];
sys_ref_r_org = [sys_ref_r_org {sys_red}];

% omega = 10.^(2:0.0005:5);
% sigma(sys_red, omega,'k');
% xlim([2*10^2 3*10^4])
% ylim([-100 20])

clear sys sys_org sys_red M_sys D_sys K_sys B_sys C_sys o E A B C V W


% p = 0.375
cd ModelsSolarPanel
load sys3_2.mat
cd ..

M_sys = sys.M;
K_sys = sys.K;
%D_sys = sys.D;
%D_sys = 10^(-8)*(M_sys+K_sys);
D_sys = 0.000001*K_sys + 0.000001*M_sys;
B_sys = sys.B;
C_sys = sys.Cp;

o = length(M_sys);

E = [K_sys, zeros(o,o); zeros(o,o), M_sys];
A = [zeros(o,o), K_sys; -K_sys, -D_sys]; 
B = [zeros(o,1); B_sys];
C = [C_sys , zeros(1,o)];

sys_org = sss(A,B,C,0,E);
[sys_red V W] = ModellReduktion(sys_org,q);
V_ref = [V_ref {V}];
W_ref = [W_ref {W}];
sys_ref_r_org = [sys_ref_r_org {sys_red}];

% omega = 10.^(2:0.0005:5);
% sigma(sys_red, omega,'k');
% xlim([2*10^2 3*10^4])
% ylim([-100 20])

clear sys sys_org sys_red M_sys D_sys K_sys B_sys C_sys o E A B C V W


% p = 0.5
cd ModelsSolarPanel
load sys3_3.mat
cd ..

M_sys = sys.M;
K_sys = sys.K;
%D_sys = sys.D;
%D_sys = 10^(-8)*(M_sys+K_sys);
D_sys = 0.000001*K_sys + 0.000001*M_sys;
B_sys = sys.B;
C_sys = sys.Cp;

o = length(M_sys);

E = [K_sys, zeros(o,o); zeros(o,o), M_sys];
A = [zeros(o,o), K_sys; -K_sys, -D_sys];
B = [zeros(o,1); B_sys];
C = [C_sys , zeros(1,o)];

sys_org = sss(A,B,C,0,E);
[sys_red V W] = ModellReduktion(sys_org,q);
V_ref = [V_ref {V}];
W_ref = [W_ref {W}];
sys_ref_r_org = [sys_ref_r_org {sys_red}];

% omega = 10.^(2:0.0005:5);
% sigma(sys_red, omega,'k');
% xlim([2*10^2 3*10^4])
% ylim([-100 20])

clear sys sys_org sys_red M_sys D_sys K_sys B_sys C_sys o E A B C V W


% p = 0.625
cd ModelsSolarPanel
load sys3_4.mat
cd ..

M_sys = sys.M;
K_sys = sys.K;
%D_sys = sys.D;
%D_sys = 10^(-8)*(M_sys+K_sys);
D_sys = 0.000001*K_sys + 0.000001*M_sys;
B_sys = sys.B;
C_sys = sys.Cp;

o = length(M_sys);

E = [K_sys, zeros(o,o); zeros(o,o), M_sys];
A = [zeros(o,o), K_sys; -K_sys, -D_sys];
B = [zeros(o,1); B_sys];
C = [C_sys , zeros(1,o)];

sys_org = sss(A,B,C,0,E);
[sys_red V W] = ModellReduktion(sys_org,q);
V_ref = [V_ref {V}];
W_ref = [W_ref {W}];
sys_ref_r_org = [sys_ref_r_org {sys_red}];

% omega = 10.^(2:0.0005:5);
% sigma(sys_red, omega,'k');
% xlim([2*10^2 3*10^4])
% ylim([-100 20])

clear sys sys_org sys_red M_sys D_sys K_sys B_sys C_sys o E A B C V W


% p = 0.75
cd ModelsSolarPanel
load sys3_5.mat
cd ..

M_sys = sys.M;
K_sys = sys.K;
%D_sys = sys.D;
%D_sys = 10^(-8)*(M_sys+K_sys);
D_sys = 0.000001*K_sys + 0.000001*M_sys;
B_sys = sys.B;
C_sys = sys.Cp;

o = length(M_sys);

E = [K_sys, zeros(o,o); zeros(o,o), M_sys];
A = [zeros(o,o), K_sys; -K_sys, -D_sys];
B = [zeros(o,1); B_sys];
C = [C_sys , zeros(1,o)];

sys_org = sss(A,B,C,0,E);
[sys_red V W] = ModellReduktion(sys_org,q);
V_ref = [V_ref {V}];
W_ref = [W_ref {W}];
sys_ref_r_org = [sys_ref_r_org {sys_red}];

% omega = 10.^(2:0.0005:5);
% sigma(sys_red, omega,'k');
% xlim([2*10^2 3*10^4])
% ylim([-100 20])

clear sys sys_org sys_red M_sys D_sys K_sys B_sys C_sys o E A B C V W
clear q