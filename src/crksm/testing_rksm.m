clear 
clc
clearvars -global

%% testdaten My_rksm

% load any benchmark system

% other models: rail_1357, rail_5177, iss, gyro, building, CDplayer, beam,
%               eady, heat-cont, LF10, random, bips98_606,
%               SpiralInductorPeec, fom

%%
sys = sss('fom')
sys = sss('rail_1357');
Opts.method = 'hammarling'
[S_P,R_Q] = lyapchol(sys);
[R_Q_trans] = lyapchol(sys.');
norm(R_Q - R_Q_trans)

Opts.method = 'adi'
[S_P,R_Q] = lyapchol(sys,Opts);
[R_Q_trans] = lyapchol(sys.',Opts);
norm(R_Q - R_Q_trans)

Opts.method = 'crksm';
% Opts.initShiftsStrategy = 'ADI'
Opts.nShifts = 10;
Opts.initShiftsStrategy = 'eigs'
% Opts.shifts = 'fixedCyclic'
% Opts.restolLyap = 1e-2;
% Opts.initShiftsStrategy = 'constant'
[S_P,R_Q,dataLyap] = lyapchol(sys,Opts);
% [R_Q_trans] = lyapchol(sys.',Opts);
% norm(R_Q - R_Q_trans)

clear
% Opts.initShiftsStrategy = 'ADI';
Opts.initShiftsStrategy = 'eigs'
Opts.nShifts = 10;
Opts.maxiter = 35;
% Opts.shifts = 'fixedCyclic'
% Opts.restolLyap = 1e-2;
% s0_inp = initializeShifts(sys,Opts.nShifts,1,Opts); %s0_inp = double(single(s0_inp));
[s0_inp,~,s0_out] = initializeShifts(sys,Opts.nShifts,1,Opts);
[sysrCrksm,V,W,S,dataCrksm] = crksm(sys,[],s0_out,Opts);

% s0_inp = initializeShifts(sys.',Opts.nShifts,1,Opts); %s0_inp = double(single(s0_inp));
% [sysrCrksmDual,VDual,WDual,SDual,dataDua] = crksm(sys.',s0_inp,Opts);

%%
clear 
sys = sss('fom');
s0_inp = [0 Inf];
Opts.purpose = 'MOR';
[sysrCrksm,V,W,S,data] = crksm(sys, s0_inp, Opts);
figure; bode(sys,'k',sysrCrksm,'r--');

clear 
sys = sss('fom');
s0_inp = [0 Inf];
[sysrCrksm,V,W,S,data] = crksm(sys, s0_inp);
figure; bode(sys,'k',sysrCrksm,'r--');

clear
sys = sss('eady');
s0_inp = [0 Inf];
Opts.purpose = 'MOR'; 
Opts.restolMOR =  1e-3;
[sysrCrksm,V,W,S,data] = crksm(sys, s0_inp, Opts);
figure; bode(sys,'k',sysrCrksm,'r--');

clear
sys = sss('iss');
s0_inp = [0 Inf];
s0_out = [0 10];
Opts.purpose = 'MOR'; 
Opts.restolMOR =  1e-3;
Opts.maxiter = 20;
[sysrCrksm,V,W,S,data] = crksm(sys, s0_inp, s0_inp, Opts);
figure; bode(sys,'k',sysrCrksm,'r--');

clear
sys = sss('rail_1357');
s0_inp = [0 Inf];
Opts.purpose = 'MOR'; 
Opts.restolMOR =  1e-3;
Opts.maxiter = 20;
[sysrCrksm,V,W,S,data] = crksm(sys, s0_inp, Opts);
figure; bode(sys,'k',sysrCrksm,'r--');

clear
sys = sss('rail_1357');
s0_inp = [0 0]; Opts.shifts = 'fixedCyclic';
Opts.purpose = 'MOR'; 
Opts.restolMOR =  1e-3;
Opts.maxiter = 20;
[sysrCrksm,V,W,S,data] = crksm(sys, s0_inp, Opts);
figure; bode(sys,'k',sysrCrksm,'r--');

%% fom model, Lyapunov
sys = loadSss('fom');

Opts.initShiftsStrategy = 'ADI';
[s0_inp,~,~,~] = initializeShifts(sys,10,1,Opts);

Opts.maxiter = 200;
%Opts.shifts = 'fixedCyclic';
Opts.restolLyap   =  1e-8; 
Opts.stopCrit  = 'residualLyap'; 
Opts.shifts = 'dynamical';
Opts.getShiftsStrategy = 'adaptive';
%Opts.lowrank = 1;
Opts.infoLyap = 1;
%Opts.rksmnorm = 'fro';
%Opts.stopCrit = 'sysr';
%Opts.equation = 'observe';
%Opts.stopCrit  = 'normChol';
%Opts.purpose = 'MOR';
[Scrksm_fom1,R,data] = lyapchol(sys,Opts);
[Scrksm_fom1,data] = lyapchol(sys',Opts);


tic
Opts.method = 'crksm';
%[Scrksm_fom1,data] = lyapchol(sys,Opts);
%[Scrksm_fom1,Rcrksm_fom1,datacrksm] = lyapchol(sys,Opts);
%[Scrksm_fom2,Rcrksm_fom2] = lyapchol(sys,Opts);
[sysr_fom,Vcrksm_fom,Wcrksm_fom,Scrksm_fom,data_crksm_fom] = crksm(sys,s0_inp,Opts);
%[sysr_fom,Vcrksm_fom,Wcrksm_fom,Rcrksm_fom,data_crksm_fom] = crksm(sys,[],s0_inp,Opts);
time_crksm_fom=toc;
clear s0_inp s0_out

% mess-adi options
Opts.method = 'adi';
Opts.messPara = 'projection';    % only for MESS
Opts.rctol = 0;
Opts.restol = 1e-8;
Opts.norm = 2;

tic
[Sadi_fom,r,data] = lyapchol(sys,Opts);
time_adi_fom = toc;

% perform tbr
tic
%sysr_crksm = tbr(sys,10,Scrksm_fom,Rcrksm_fom);
time_tbr_crksm = toc;

tic
%sysr_adi = tbr(sys,10,Sadi_fom,Radi_fom);
time_tbr_adi = toc;
% call rk
%[sysr_rk,Vrk,Wrk] = rk(sys,s0_inp);
% test subspace
%alpha1 = subspace(Vcrksm,Vrk);
clear Opts

%% CDplayer model, MOR
sys = loadSss('CDplayer');

Opts.strategy = 'ADI';
[s0_inp,~,~,~] = initializeShifts(sys,10,1,Opts);

Opts.maxiter = 200;
%Opts.shifts = 'fixedCyclic';
%Opts.restolLyap   =  1e-8; 
%Opts.stopCrit  = 'residualLyap'; 
Opts.shifts = 'dynamical';
Opts.strategy = 'adaptive';
%Opts.rksmnorm = 'fro';
%Opts.stopCrit = 'sysr';
%Opts.equation = 'observe';
%Opts.stopCrit  = 'normChol';
Opts.purpose = 'MOR';

tic
%[Scrksm,Rcrksm] = lyapchol(sys,Opts);
%[sysr_CD,Vcrksm_CD,Wcrksm_CD,Scrksm_CD,data_crksm_CD] = crksm(sys,s0_inp,Opts);
time_crksm_CD=toc;
clear s0_inp s0_out

% initialize irka
% s0_inp = zeros(1,10);
% Rt = ones(m,size(s0_inp,10));
% Lt = ones(p,size(s0_inp,10));

% tic
% [sysr_irka_CD, ~, ~, s0_inp, Rt, Lt] = irka(sys,s0_inp,Rt,Lt);
% time_irka_CD = toc;

% call rk
%[sysr_rk,Vrk,Wrk] = rk(sys,s0_inp);
% test subspace
%alpha1 = subspace(Vcrksm,Vrk);
clear Opts

%% Rail 1357 Lyapunov block
sys = loadSss('rail_1357');

Opts.strategy = 'ADI';
[s0_inp,~,~,~] = initializeShifts(sys,10,1,Opts);

Opts.maxiter = 200;
%Opts.shifts = 'fixedCyclic';
Opts.restolLyap   =  1e-8; 
Opts.stopCrit  = 'residualLyap'; 
Opts.shifts = 'dynamical';
Opts.strategy = 'adaptive';
%Opts.rksmnorm = 'fro';
%Opts.stopCrit = 'sysr';
%Opts.equation = 'observe';
%Opts.stopCrit  = 'normChol';


tic
%[Scrksm,Rcrksm] = lyapchol(sys,Opts);
%[sysr_rail13,Vcrksm_rail13,Wcrksm_rail13,Scrksm_rail13,data_crksm_rail13] = crksm(sys,s0_inp,Opts);
time_crksm_rail13=toc;
clear s0_inp s0_out

% mess-adi options
Opts.method = 'adi';
Opts.messPara = 'projection';    % only for MESS
Opts.rctol = 0;
Opts.restol = 1e-8;
Opts.norm = 2;

tic
%[Sadi_rail13] = lyapchol(sys,Opts);
time_adi_rail13 = toc;
clear Opts

%% Rail 5177 Lyapunov block
sys = loadSss('rail_5177');

Opts.strategy = 'ADI';
[s0_inp,~,~,~] = initializeShifts(sys,10,1,Opts);

Opts.maxiter = 200;
%Opts.shifts = 'fixedCyclic';
Opts.restolLyap   =  1e-8; 
Opts.stopCrit  = 'residualLyap'; 
Opts.shifts = 'dynamical';
Opts.strategy = 'adaptive';
%Opts.rksmnorm = 'fro';
%Opts.stopCrit = 'sysr';
%Opts.equation = 'observe';
%Opts.stopCrit  = 'normChol';


tic
%[Scrksm,Rcrksm] = lyapchol(sys,Opts);
%[sysr_rail51,Vcrksm_rail51,Wcrksm_rail51,Scrksm_rail51,data_crksm_rail51] = crksm(sys,s0_inp,Opts);
time_crksm_rail51=toc;
clear s0_inp s0_out

% mess-adi options
Opts.method = 'adi';
Opts.messPara = 'projection';    % only for MESS
Opts.rctol = 0;
Opts.restol = 1e-8;
Opts.norm = 2;

tic
%[Sadi_rail51] = lyapchol(sys,Opts);
time_adi_rail_51 = toc;
clear Opts

%% Rail 5177 Lyapunov tangential
sys = loadSss('rail_5177');

Opts.strategy = 'eigs';
[s0_inp,Rt,~,~] = initializeShifts(sys,10,1,Opts);

Opts.maxiter = 200;
%Opts.shifts = 'fixedCyclic';
Opts.restolLyap   =  1e-8; 
Opts.stopCrit  = 'residualLyap'; 
Opts.shifts = 'dynamical';
Opts.strategy = 'adaptive';
%Opts.rksmnorm = 'fro';
%Opts.stopCrit = 'sysr';
%Opts.equation = 'observe';
%Opts.stopCrit  = 'normChol';


tic
%[Scrksm,Rcrksm] = lyapchol(sys,Opts);
[sysr_rail51t,Vcrksm_rail51t,Wcrksm_rail51t,Scrksm_rail51t,data_crksm_rail51t] = crksm(sys,s0_inp,Rt,Opts);
time_crksm_rail51t=toc;
clear s0_inp s0_out Opts

%% gyro noch machen


%% testing quality of solution (when dimension is high, comment out n>1500)

% with low rank factor
P = Scrksm*Scrksm';
YS = sys.A*P*sys.E' + sys.E*P*sys.A' + sys.B*sys.B';
YS_norm = norm(YS);

Q = Rcrksm*Rcrksm';
YR = sys.A'*Q*sys.E + sys.E'*Q*sys.A + sys.C'*sys.C;
YR_norm = norm(YR);

% no low rank factor
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
     P = Vcrksm*Pr*Vcrksm';
     YS = A*P*E' + E*P*A' + B*B';
     YS_norm = norm(YS);
     Qr = R'*R;
     Q = Wcrksm*Qr*Wcrksm';
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


