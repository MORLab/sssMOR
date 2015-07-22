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
load gyro; %build, beam, fom, rail_1357, rail_1357
if ~exist('A','var') && exist('M','var') %2nd order
    E = blkdiag(M,M);
    A = [zeros(size(M)),M; - K, -1e-6*K];
    B = [zeros(size(B)); B];
    C = [C, zeros(size(C))];
    clear M K
end
    
if size(B,2)>1, B = B(:,1);end
if size(C,1)>1, C = C(1,:);end
if ~exist('D','var'), D = zeros(size(C,1),size(B,2)); end
if ~exist('E','var'), E = speye(size(A)); end
    
sys = sss(A,B,C(1,:),D,E);

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

if size(A,1) < 2e3
    figure;bode(ss(sys),'b-',ss(sysr),'r--',ss(sysrRitz),'g-.');
    h2Norm      = norm(sys-sysr)
    h2NormRitz  = norm(sys-sysrRitz)
else
    figure;bode(sysr,'r--');
    hold on; bode(sysrRitz,'g-.');
    h2Norm = norm(sysr)
    h2NormRitz = norm(sysrRitz)
end

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

