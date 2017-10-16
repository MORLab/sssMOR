clear all
clc
clearvars -global

%% testdaten My_rksm

% load any benchmark system
sys = loadSss('fom');

% other models: rail_1357, rail_5177, iss, gyro, building, CDplayer, beam,
%               eady, heat-cont, LF10, random, bips98_606,
%               SpiralInductorPeec, fom

% declare system matices
A = sys.A;
B = sys.B;
C = sys.C;
E = sys.E;
D = sys.D;

sys_new = sss(A,B,C,D,E);
m = size(B,2);
p = size(C,1);



%% testing function for several cases

% initialize shifts with irka or function 'initializeShifts'
% irka
s0_inp = zeros(1,10);
Rt = ones(m,size(s0_inp,10));
Lt = ones(p,size(s0_inp,10));
tic
[sysr_irka, ~, ~, s0_inp, Rt, Lt] = irka(sys,s0_inp,Rt,Lt);
time_irka = toc;

Opts.strategy = 'eigs';
[s0_inp,Rt,s0_out,Lt] = initializeShifts(sys,10,1,Opts);
%s0_inp = reshape(s0,[1,20]);
%Opts.orth = 0;
%[V_arn,Sv_arn,Rv_arn] = arnoldi(sys.E,sys.A,sys.B,s0_inp,Opts);
%s0_inp = [0 inf 0 inf 0 inf 0 inf 0 inf 0 inf];
%s0_out = s0;


% call crksm and other functions
Opts.maxiter_rksm = 40;
Opts.rctolLyap    =  1e-7; 
Opts.residual  = 'residual_lyap';  
%Opts.residual  = 'norm_chol';
%Opts.purpose = 'MOR';
Opts.shifts = 'dynamical';
%Opts.shifts = 'fixedCyclic';
Opts.strategy = 'adaptive';
%Opts.rksmnorm = 'fro';
Opts.stopCrit = 'sysr';
%Opts.equation = 'observe';

    % call rk
    tic
    %[sysr_rk,V_rk,W_rk] = rk(sys,s0_inp);
    time_rk = toc;

    % call function crksm
    tic
    [sysr,V_crksm,W_crksm,S,R,data_crksm] = crksm(sys,s0_inp,Opts);
    time_crksm=toc;
      
    % test for tangential cascaded subspace
    %tic
    %[V_tan_cas,W_tan_cas] = test1(sys,s0_inp,Rt,Opts);
    %time_test1 = toc;
        
    alpha1 = subspace(data_out.out4,V_rk);
    alpha2 = subspace(data_out.out4,Varnoldi);
    alpha3 = subspace(Varnoldi,V_rk);



%% testing quality of solution (when dimension is high, comment out n>1500)
% V = data_out.V_basis;
% %W = data_out.W_basis;
% if isempty(W)
%     Pr = S'*S;
%     P = V*Pr*V';
%    % P  = Pr'*Pr;
% else
%     Pr = W*S*V';
%     P  = Pr*Pr';
% end
%  
% Y = A*P*E' + E*P*A' + B*B';
% Y_norm = norm(Y);


% comparing with other hammarling-method
% Opts.method = 'hammarling';
% [Shammarling,Rhammarling] = lyapchol(sys,Opts);
% Phammarling = Shammarling*Shammarling';
% Yhammarling = A*Phammarling*E' + E*Phammarling*A' + B*B';
% Yhammarling_norm = norm(Yhammarling);



% mess-adi options
% Opts.method = 'adi';
% Opts.messPara = 'heur';    % only for MESS
% Opts.rctol = 0;
% Opts.restol = 1e-5;
% Opts.norm = 2;
% 
% call 
% tic
% [Sadi,Radi] = lyapchol(sys,Opts);
% time_adi = toc;
% 
% testing quality of solution (when dimension is high, comment out n>1500)
% Padi = Sadi*Sadi';
% Yadi = A*Padi*E' + E*Padi*A' + B*B';
% Yadi_norm = norm(Yadi);

% save data from workspace in your path, comment the following in
% file_name = 'eady';
% save(file_name)

% genrate plot (if you want a plot of the residual norm comment the following in ;))

% plot(1:size(S,2)/m,data_out.norm_val, '-ob','LineWidth',2);
% hold on
% plot(1:size(Sadi,2)/m,Radi.res,'-+r','LineWidth',2);
% title('Convergence Rate "fom" model CRKSM ');
% xlabel('Number of Iterations');
% ylabel('Convergence Rate');
% legend('CRKSM','Location','southwest');
% legend('CRKSM','ADI','Location','southwest');


% generate a filename.tikz file for using matlab2tikz in latex,
% matlab2tikz-directory must be in yout path (I removed this directory, so you have to add it again) 
% note: cleanfigure-command is helpful!

% cleanfigure
% matlab2tikz( 'filename.tikz' );


