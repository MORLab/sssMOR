clear all
clc
clearvars -global

%% testdaten My_rksm

% load any benchmark system
sys = loadSss('rail_1357');

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
s0_inp = zeros(1,20);
Rt = ones(m,size(s0_inp,20));
Lt = ones(p,size(s0_inp,20));
tic
%[sysr_irka, ~, ~, s0_inp, Rt, Lt] = irka(sys,s0_inp,Rt,Lt);
time_irka = toc;

Opts.strategy = 'ROM';
[s0_inp1,Rt1,s0_out1,Lt1] = initializeShifts(sys,10,1,Opts);
%s0_inp1 = initializeShifts(sys,10,1,Opts);
Opts.strategy = 'ADI';
[s0_inp2,Rt2,s0_out2,Lt2] = initializeShifts(sys,10,1,Opts);
%s0_inp2 = initializeShifts(sys,10,1,Opts);
%s0_inp2 = s0_inp2';
Opts.strategy = 'ROM';




% call crksm and other functions
Opts.maxiter = 200;
Opts.restolLyap   =  1e-6; 
Opts.stopCrit  = 'residualLyap';  
%Opts.stopCrit  = 'normChol';
%Opts.purpose = 'MOR';
Opts.shifts = 'dynamical';
%Opts.shifts = 'fixedCyclic';
Opts.strategy = 'adaptive';
Opts.adiMethod = 'heur';
%Opts.rksmnorm = 'fro';
%Opts.stopCrit = 'sysr';
%Opts.equation = 'observe';

    % call rk
    tic
    %[sysr_rk,Vrk,Wrk] = rk(sys,s0_inp);
    time_rk = toc;

    % call function crksm
    tic
    [sysr,Vcrksm,Wcrksm,S,R,data_crksm] = crksm(sys,[],s0_inp2,Opts);
    time_crksm=toc;
    
    % mess-adi options
    Opts.method = 'adi';
    Opts.messPara = 'heur';    % only for MESS
    Opts.rctol = 0;
    Opts.restol = 1e-6;
    Opts.norm = 2;
    Opts.info = 1;
    
    tic
    [Sadi,Radi] = lyapchol(sys,Opts);
    time_adi = toc;
    
    % test subspace
    %alpha1 = subspace(Vcrksm,Vrk);

%% testing quality of solution (when dimension is high, comment out n>1500)

if isempty(R)
    Pr = S'*S;
    P = Vcrksm*Pr*Vcrksm';
    YS = A*P*E' + E*P*A' + B*B';
    YS_norm = norm(YS);
elseif isempty(S)
    Qr = R'*R;
    Q = Wcrksm*Qr*Wcrksm';
    YR = A'*Q*E + E'*Q*A + C'*C;
    YR_norm = norm(YR);
else
     Pr = S'*S;
     P = Wcrksm*Pr*Vcrksm';
     YS = A*P*E' + E*P*A' + B*B';
     YS_norm = norm(YS);
     Qr = R'*R;
     Q = Wcrksm*Qr*Vcrksm';
     YR = A'*Q*E + E'*Q*A + C'*C;
     YR_norm = norm(YR);
end
    



% ADI solution
Padi = Sadi*Sadi';
YSadi = A*Padi*E' + E*Padi*A' + B*B';
YSadi_norm = norm(YSadi);

Qadi = Radi*Radi';
YRadi = A'*Qadi*E + E'*Qadi*A + C'*C;
YRadi_norm = norm(YRadi);

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


