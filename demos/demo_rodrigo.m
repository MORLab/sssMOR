%demo rodrigo
clear, clc, close all
load build
% load iss

sys = sss(A,B,C,0);

%% Bode
bode(sys); hold on
% [mag2,ph2,w2] = bode(ss(sys));

%% Reduktion
n_red = 10;
%Art 1
s0 = zeros(1,n_red);
%Art 2
% s0 = [0 1;n_red/2, n_red/2];

% sysr1 = RK(sys,s0);
% bode(sysr1,'--g')

sysr2 = RK(sys,[],s0);
bode(sysr2,'-.r')

sysr3 = RK(sys,s0,s0);
bode(sysr3,'-k')


sysr4 = IRKA(sys,s0);
bode(sysr4,'-.c')

opts = {'complete',0};
sysr5 = IRKA_analyze(sys, s0,[],[],opts);
bode(sysr5,'-.b')


