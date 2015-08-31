%% Define version/functions to test

% dirname = spritf(...

% go to target directory
% cd()

%% Test Bode
%1) bode plot
sys = loadSss('build');

bode(sys)
%% Test RK and Arnoldi
clear, 
clc

load build
sys = sss(A,B,C);

% s0 = rand(1,20);
s0 = rand(1,10) + 1i*rand(1,10);
s0 = [s0, conj(s0)];
tic,[sysr1i,V1i,W1i]    = rk(sys,s0);t1i = toc;
tic,[sysr1o,V1o,W1o]    = rk(sys,[],s0); t1o = toc;
tic,[sysr2,V2,W2]       = rk(sys,s0,s0); t2 = toc;

norm(V1i-V2)
norm(W1o-W2)

sysr3 = sss(W1o'*sys.a*V1i,W1o'*sys.b, sys.c*V1i,sys.d,W1o'*sys.e*V1i);
figure;bode(sys,'k');hold on, bode(sysr2,'r-'); bode(sysr3,'g--');

norm(cplxpair(eig(sysr2))-cplxpair(eig(sysr3)))
norm(sysr3-sysr2)/norm(sysr3)
%%  Test Arnoldi (invsolve)
clear, clc
n = 1e2; density = 0.1;
A = sprand(n,n,density)+ 1i*sprand(n,n,density); 
E = sprand(n,n,density); b = sprand(n,1,density) + 1i*sprand(n,1,density);
Id = speye(size(A));

[L,U,p,o,S]=lu(A,'vector');
[Lt,Ut,pt,ot,St] = lu(A.','vector');
    
% Input Krylov
fprintf('--- input Krylov ---');
    %1st direction
    xCorr = A\b;
    %a) invsolve
    tic
    xInvsolve=S\b;
    xInvsolve=L\xInvsolve(p,:);
    xInvsolve(o,:)=U\xInvsolve;
    tInvsolve = toc
    eInvsolve = norm(xCorr-xInvsolve)
    
    %b) direct computation
    tic
    xDirect = Id(:,o)*(U\(L\(Id(:,p)'*(S\b))));
    tDirect = toc
    eDirect = norm(xCorr-xDirect)
    
    %c) New: a mixture of both a and b
    tic
    xNew(o,:) = U\(L\(S(:,p)\b));
    tNew = toc
    eNew = norm(xCorr-xNew)

    
    %2nd direction
    xCorr = A\(E*xCorr);
    %a) invsolve
    xInvsolve=S\(E*xInvsolve);
    xInvsolve=L\xInvsolve(p,:);
    xInvsolve(o,:)=U\xInvsolve;
    eInvsolve = norm(xCorr-xInvsolve)
    
    %b) direct computation
    xDirect = Id(:,o)*(U\(L\(Id(:,p)'*(S\(E*xDirect)))));
    eDirect = norm(xCorr-xDirect)
    
    %c) New: a mixture of both a and b
    xNew(o,:) = U\(L\(S(:,p)\(E*xNew)));
    eNew = norm(xCorr-xNew)
    
% Output Krylov
fprintf('--- output Krylov ---');
    %1st direction
    xCorr = A.'\b;
    xCorr2(ot,:) = Ut\(Lt\(St(:,pt)\b));
    eLUt = norm(xCorr-xCorr2)

    %a) direct computation
    tic
    xDirect = (S*Id(:,p)).'\(L.'\(U.'\(Id(:,o)'*b)));
    tDirect = toc
    eDirect = norm(xCorr-xDirect)
    
    %c) New: a mixture of both a and b
    tic
    xNew = (S(:,p)).'\(L.'\(U.'\(b(o,:))));
    tNew = toc
    eNew = norm(xCorr-xNew)
    
    %2nd direction
    xCorr2(ot,:) = Ut\(Lt\(St(:,pt)\(E'*xCorr)));
    xCorr = A.'\(E'*xCorr);
    eLUt = norm(xCorr-xCorr2)
    
    %b) direct computation
    xDirect = (S*Id(:,p)).'\(L.'\(U.'\(Id(:,o)'*(E'*xDirect))));
    eDirect = norm(xCorr-xDirect)
    
    %c) New: a mixture of both a and b
    temp = E'*xNew;
    xNew = (S(:,p)).'\(L.'\(U.'\(temp(o,:))));
    eNew = norm(xCorr-xNew)  
%%  Test IRKA 
%   [Ale]
clear, clc

load build
sys = sss(A,B,C);
s0 = 100*rand(1,40);
Opts = struct('maxiter',100,'epsilon',1e-3,'stopCrit','combAll','verb',1);
[sysr, V, W, s0, s0_traj] = irka(sys, s0, Opts);
analyzeMor(sys,sysr);
%%  Test CURE
close all, clear, clc

% NOTE: obsolete --v
% load beam; %build, beam, fom, rail_1357, rail_1357, gyro
% if ~exist('A','var') && exist('M','var') %2nd order
%     E = blkdiag(M,M);
%     A = [zeros(size(M)),M; - K, -1e-6*K];
%     B = [zeros(size(B)); B];
%     C = [C, zeros(size(C))];
%     clear M K
% end
%     
% if size(B,2)>1, B = B(:,1);end
% if size(C,1)>1, C = C(1,:);end
% if ~exist('D','var'), D = zeros(size(C,1),size(B,2)); end
% if ~exist('E','var'), E = speye(size(A)); end
% 
% sys = sss(A,B,C(1,:),D,E);
% -- ^

sys = loadSss('beam');
   
Opts.cure.init = 'slm';
Opts.cure.test = 0;
Opts.spark.test = 0;
Opts.cure.stopval = 20;

Opts.verbose = 0; %show progress text?
% Opts.cure.red = 'rk'; %reduction algorithm
% Opts.cure.nk = 10; % reduced order at each step

Opts.mespark.ritz = 0;
tic, sysr = cure(sys,Opts); tCure = toc

Opts.mespark.ritz = 1;
% Opts.spark.test = 1;
tic, sysrRitz = cure(sys,Opts); tCureRitz = toc 

if size(sys.A,1) < 2e3
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

%%  Test all MOR functions
clc
fprintf('------------------------------------\n');
fprintf('Test all MOR functions:\n');
fprintf('------------------------------------\n');

benchmarkName = {'build','beam'}; %to be expanded
morFunction = {'rk','tbr','modalMor','irka','cure'};

for iMorFun = 1:length(morFunction)
    % define a handle for the reduction
    switch morFunction{iMorFun}
        case 'rk'
            q = 10;
            s0 = zeros(1,q);
            redFun = @(sys) rk(sys,s0,s0);
        case 'tbr'
            q = 10;
            redFun = @(sys) tbr(sys,q);
        case 'modalMor'
            q = 10;
            redFun = @(sys) modalMor(sys, q);
        case 'irka'
            q = 10;
            redFun = @(sys) irka(sys,zeros(1,q));
        case 'cure'
            redFun = @(sys) cure(sys);
        otherwise
            error('This reduction method is not defined for testing');
    end
    for jBenchmark = 1:length(benchmarkName)
        fprintf('Testing reduction of %s with %s...\n',...
            benchmarkName{jBenchmark}, morFunction{iMorFun});
        
        %TODO: if several reductions should be conducted for one reduction
        %methods, pass a CELL ARRAY of functions handles and go through all
        %of them!
        
        sys  = loadSss(benchmarkName{jBenchmark});
        sysr = redFun(sys); 
        analyzeMor(sys,sysr,struct('mode','red'));
        pause(3), close
    end
end
fprintf('------------------------------------\n');
fprintf('Test all MOR functions finished!\n')
fprintf('------------------------------------\n');