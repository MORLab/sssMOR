%% Define version/functions to test

dirname = spritf(...

% go to target directory
cd()

%% Test Bode
%1) bode plot
load(build)
sys = sss(A,B,C,D,E)

bode(sys)

%% Test CURE
close all, clear, clc
load fom; %beam, fom
sys = sss(A,B,C);

Opts.CURE.init = 'slm';
Opts.CURE.test = 0;
Opts.SPARK.test = 0;
Opts.CURE.stopval = 20;

Opts.verbose = 0; %show progress text?
% Opts.CURE.red = 'RK'; %reduction algorithm
% Opts.CURE.nk = 10; % reduced order at each step

Opts.MESPARK.ritz = 0;
tic, sysr = CURE(sys,Opts); tCure = toc

Opts.MESPARK.ritz = 1;
% Opts.SPARK.test = 1;
tic, sysrRitz = CURE(sys,Opts); tCureRitz = toc 

figure;bode(ss(sys),'b-',ss(sysr),'r--',ss(sysrRitz),'g-.');
h2Norm      = norm(sys-sysr)
h2NormRitz  = norm(sys-sysrRitz)

%  Test initializeShifts (within CURE)
% for ii = 1:20
%     Opts.CURE.nk = ii;
%     Opts.CURE.init = 1;
%     Opts.zeroThres = 0;
%     s0 = initializeShifts(sys,Opts);
%     plot(real(s0),imag(s0),'*')
%     title(sprintf('Opts.CURE.nk = %i',ii));
%     pause
% end

