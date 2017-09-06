%% Testing Cascaded Tangential Krylov
%
clear
clc


% get a random system matrices

Opts.orth = '2mgs';     

[rk_stand,rk_cascad,Vstand,Vcascad,alpha,aa,bb] = test1(Opts);

sys = loadSss('iss');
A = sys.A;
B = sys.B;
C = sys.C;
E = sys.E;

m = size(B,2);      n = size(A,2);
p = size(C,1);


% compute shifts using irka
shifts = zeros(1,10);
Rt = ones(m,size(shifts,10));
Lt = ones(p,size(shifts,10));
[~, ~, ~, s0, T, Lt] = irka(sys, shifts,Rt,Lt);
%T = ones(m,size(s,2));

[sysr, Vrk, Wrk] = rk(sys,s0);


%W = data_out.W_basis;
alpha=subspace(Vcascad,Vrk);
beta=subspace(Vstand,Vrk);

