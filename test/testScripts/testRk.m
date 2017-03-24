classdef testRk < sssTest
% testRk - testing of rk.m
%
% Description:
%   The function rk.m is tested (5 tests) on:
%    + comparing sysr to a directly calculated solution based on arnoldi.m
%    + one-sided (V=W), double-sided 
%    + Ar, Er purely real
%    + rank(Ar), rank(Er) full
%    + Neither Inf nor NaN in Ar, Er
%    + s0: purely real, purely imaginary, zero, Inf 
%    + test systems: building, beam, random, SpiralInductorPeec 
%      (with E-matrix), LF10 (with E-matrix).
%    + is moment matching achieved
%    + are the matrices of Sylvester equation correct (B_, Rv, C_, Lw)
%
%    + TO DO: Error in testRk5
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
%               Lisa Jeschek
% Last Change:  06 Apr 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------
    
    methods(Test)
         function testRk1 (testCase) 
              %one-sided reduction isempty(s0_out), s0: real (multiple value)
              load('building.mat'); sys = sss(A,B,C);
              n = 10; s0val = 100; s0 = ones(1,n)*s0val; 
              
              [sysr, V, W, B_, ~, Rv, C_, ~, Lw] = rk(sys, s0);
              actSolution={full(sysr.A), full(sysr.B), full(sysr.C), V, W, ...
                  Rv, B_};
              
              expV = arnoldi(speye(size(A)),A,B, s0);
              [expRv,expB_] = getSylvester(sys,sysr,V);
              expSolution={expV'*A*expV, expV'*B, C*expV, expV, expV,...
                   expRv, expB_}; 
               
              % check for moment matching as well
              actM = moments(sysr,s0val,n); actSolution = [actSolution,{actM}];
              expM = moments(sys,s0val,n); expSolution = [expSolution,{expM}];
               
              verification(testCase, actSolution, expSolution, sysr);
              verifyEmpty(testCase, C_, ...
                    'C_ is not empty');
              verifyEmpty(testCase, Lw,...
                    'Lw is not empty');
              verifyEqual(testCase, V, W,...
                    'V is not equal W');
         end
         function testRk2 (testCase) 
              %one-sided reduction isempty(s0_in), s0: imag (multiple value)
              load('beam.mat'); sys = sss(A,B,C);
              n = 5; s0val = 100; s0 = [ones(1,n)*s0val*1i,-ones(1,n)*s0val*1i]; 
              
              [sysr, V, W, B_, ~, Rv, C_, ~, Lw] = rk(sys, [], s0);
              actSolution={full(sysr.A), full(sysr.B), full(sysr.C), V, W, Lw, C_};
              
              expW = arnoldi(speye(size(A)),A',C',s0);
              [expLw,expC_] = getSylvester(sys,sysr,W,'W');
              expSolution={expW'*A*expW, expW'*B, C*expW, expW, expW, expLw,...
                    expC_};
                
              % check for moment matching as well
              actM = moments(sysr,[s0(1),s0(end)], [n,n]); actSolution = [actSolution,{actM}];
              expM = moments(sys,[s0(1),s0(end)], [n,n]); expSolution = [expSolution,{expM}];
                
              verification(testCase, actSolution, expSolution, sysr);
              verifyEmpty(testCase, B_, ...
                    'B_ is not empty');
              verifyEmpty(testCase, Rv,...
                    'Rv is not empty');
              verifyEqual(testCase, V, W,...
                    'V is not equal W');
         end         
         function testRk3 (testCase) 
              %two-sided reduction without E-matrix, s0: zero, imag, real, Inf
              load('random.mat'); sys = sss(A,B,C); IP = @(x,y) x'*y;
              s0 = [1+4i, 1-4i,4+4i, 4-4i, 0, 0, 5, 300, Inf,  Inf];
              s0moment = s0([1,3,5,7:9]); n = [2, 2, 4, 2, 2, 4];
              
              % careful: Sylvester EQ probably does not hold for shifts at
              %         Infinity!!
              [sysr, V, W] = rk(sys,s0 , s0, IP);
              actSolution={full(sysr.A), full(sysr.B), full(sysr.C), V, W};
              
              [expV,~,~,expW,~] = arnoldi(speye(size(A)),A,B, C, s0, IP);
              expSolution={expW'*A*expV, expW'*B, C*expV, expV, expW};
              
              % check for moment matching as well
              actM = moments(sysr,s0moment, n); actSolution = [actSolution,{actM}];
              expM = moments(sys,s0moment, n); expSolution = [expSolution,{expM}];
              
              % check for orthogonality
              actOrthoV = IP(V,V); actOrthoW = IP(W,W);
              verifyEqual(testCase, {actOrthoV, actOrthoW}, ...
                  {eye(size(V,2)), eye(size(W,2))},'AbsTol', 1e-6, ...
                  'Projection matrices not orthogonal');

              verification(testCase, actSolution, expSolution, sysr);
              
         end 
         function SISOTwoSidedDescriptor (testCase) 
              %two-sided reduction with E-matrix SISO
              sys = loadSss('SpiralInductorPeec');
                
                Opts.stopCrit = 's0';
                %  get good shifts
                n = 6; r = ones(sys.m,n); l = ones(sys.p,n);
                sysrIrka = irka(sys, zeros(1,n),r, l,Opts);
                s0 = -eig(sysrIrka).'; s0moment = s0; n = 2;
            
              [sysr, V, W, B_, ~, Rv, C_, ~, Lw] = rk(sys,s0,s0);              
              [expV,~,~,expW,~] = arnoldi(sys.E,sys.A,sys.B,sys.C,s0);
              
              % The transpose LU problem can be ill conditioned, check the
              % subspaces instead of the actual matrices!
              actSolution={rank([V,expV]), rank([W,expW])};
              expSolution={size(V,2), size(W,2)};
              
              % Add Sylvester EQ matrices
              [expRv,expB_] = getSylvester(sys,sysr,V);
              [expLw,expC_] = getSylvester(sys,sysr,W,'W');
              
              actSolution = [actSolution, {Rv, B_, Lw, C_}];
              expSolution = [expSolution, {expRv,expB_,expLw,expC_}];
             
              % res = norm(sys.A*V - sys.E*V*(sysr.E\sysr.A)-B_*Rv);
              % res = norm(sys.A*V - sys.E*V*(sysr.E\sysr.A)-expB_*expRv);
              
              % check for moment matching as well
              actM = moments(sysr,s0moment, n); actSolution = [actSolution,{actM}];
              expM = moments(sys,s0moment, n); expSolution = [expSolution,{expM}];
              
              verification(testCase, actSolution, expSolution, sysr);
         end 
         function MIMOoneSidedRealShifts (testCase)
             sys = loadSss('CDplayer');
             
             n = 4; s0 = 100*rand(1,n);
             
             % Input
             Rt = 100*rand(sys.m,n);
             
             sysr = rk(sys,s0,Rt);
             
             % verify moment matching
             Me = moments(sys-sysr,s0,1);
             Rt3D  = permute(Rt,[1,3,2]);
             MeR = mmat(Me,Rt3D);
             
             verifyLessThanOrEqual(testCase,abs(MeR),1e-8,'MeR: Moments do not match')   
             
             % Output
             Lt = 100*rand(sys.m,n);
             
             sysr = rk(sys,[],s0,[],Lt);
             
             % verify moment matching
             Me = moments(sys-sysr,s0,1);
             Lt3D  = permute(Lt.',[3,2,1]);
             MeL = mmat(Lt3D,Me);
             
             verifyLessThanOrEqual(testCase,abs(MeL),1e-8,'MeL: Moments do not match') 
         end
         function MIMOoneSidedComplexShiftsRealModel (testCase)
             sys = loadSss('CDplayer');
             
             n = 4; 
             s0 = 100*(rand(1,n)+1i*randn(1,n)); s0([2,4]) = conj(s0([1,3]));

             % Input
             Rt = 100*(rand(sys.m,n)+1i*randn(sys.m,n)); Rt(:,[2,4]) = conj(Rt(:,[1,3]));
             
             sysr = rk(sys,s0,Rt);
             
             % verify moment matching
             Me = moments(sys-sysr,s0,1);
             Rt3D  = permute(Rt,[1,3,2]);
             MeR = mmat(Me,Rt3D);
             
             verifyLessThanOrEqual(testCase,abs(MeR),1e-8,'MeR: Moments do not match')   
             
             % Output
             Lt = 100*(rand(sys.m,n)+1i*randn(sys.m,n)); Lt(:,[2,4]) = conj(Lt(:,[1,3]));
             
             sysr = rk(sys,[],s0,[],Lt);
             
             % verify moment matching
             Me = moments(sys-sysr,s0,1);
             Lt3D  = permute(Lt.',[3,2,1]);
             MeL = mmat(Lt3D,Me);
             
             verifyLessThanOrEqual(testCase,abs(MeL),1e-8,'MeL: Moments do not match') 
         end
         function MIMOoneSidedComplexShiftsComplexModel (testCase)
             sys = loadSss('CDplayer');
             
             n = 4; Opts.real = false;
             
             s0 = 100*(rand(1,n)+1i*randn(1,n));
            
             % Input
             Rt = 100*(rand(sys.m,n)+1i*randn(sys.m,n));
             
             sysr = rk(sys,s0,Rt,Opts);
             
             % verify moment matching
             Me = moments(sys-sysr,s0,1);
             Rt3D  = permute(Rt,[1,3,2]);
             MeR = mmat(Me,Rt3D);
             
             % NOTE due to bad conditioning of random shifts, the error
             % tolerance should not be to tight here..
             verifyLessThanOrEqual(testCase,abs(MeR),1e-3,'MeR: Moments do not match')   
             
             % Output
             Lt = 100*(rand(sys.m,n)+1i*randn(sys.m,n));
             
             sysr = rk(sys,[],s0,[],Lt,Opts);
             
             % verify moment matching
             Me = moments(sys-sysr,s0,1);
             Lt3D  = permute(Lt.',[3,2,1]);
             MeL = mmat(Lt3D,Me);
             
             verifyLessThanOrEqual(testCase,abs(MeL),1e-8,'MeL: Moments do not match') 
         end
         function MIMOtwoSidedComplexShiftsComplexModel (testCase) 
             sys = loadSss('CDplayer');
             
             n = 4; 
             %  Complex shifts and tangential directions
             Rt = ones(sys.m,n)+1i*ones(sys.m,n); 
             Lt = ones(sys.p,n)-1i*ones(sys.p,n);
             s0 = 100*(rand(1,n)+1i*rand(1,n));
             
             Opts.real = false;
             sysr = rk(sys,s0,s0,Rt,Lt,Opts);
             
             % verify moment matching
             Me2 = moments(sys-sysr,s0,2);
             Me = Me2(:,:,1:2:end); Me2(:,:,1:2:end) = [];
             
             Lt3D  = permute(Lt.',[3,2,1]);
             Rt3D  = permute(Rt,[1,3,2]);
             
             MeL = mmat(Lt3D,Me);
             MeR = mmat(Me,Rt3D);
             MeH = mmat(mmat(Lt3D,Me),Rt3D);
             
             verifyLessThanOrEqual(testCase,abs(MeL),1e-6,'MeL: Moments do not match')
             verifyLessThanOrEqual(testCase,abs(MeR),1e-6,'MeR: Moments do not match')
             verifyLessThanOrEqual(testCase,abs(MeH),1e-6,'MeH: Moments do not match')
         end
         function MIMOtwoSidedComplexShiftsRealModel (testCase) 
             % test if sorting shifts/tangential directions messes up
             % something
             
             sys = loadSss('CDplayer');
             
             n = 4; 
             %  Complex shifts and tangential directions
             Rt = 100*(rand(sys.m,n)+1i*randn(sys.m,n)); Rt(:,[2,4]) = conj(Rt(:,[1,3]));
             Lt = 100*(rand(sys.p,n)+1i*randn(sys.p,n)); Lt(:,[2,4]) = conj(Lt(:,[1,3]));
             s0 = 100*(rand(1,n)+1i*randn(1,n)); s0([2,4]) = conj(s0([1,3]));
             
             Opts.real = true; 
             sysr = rk(sys,s0,s0,Rt,Lt,Opts);
             
             % verify moment matching
             Me2 = moments(sys-sysr,s0,2);
             Me = Me2(:,:,1:2:end); Me2(:,:,1:2:end) = [];
             
             Lt3D  = permute(Lt.',[3,2,1]);
             Rt3D  = permute(Rt,[1,3,2]);
             
             MeL = mmat(Lt3D,Me);
             MeR = mmat(Me,Rt3D);
             MeH = mmat(mmat(Lt3D,Me2),Rt3D);
             
             verifyLessThanOrEqual(testCase,abs(MeL),1e-5,'MeL: Moments do not match')
             verifyLessThanOrEqual(testCase,abs(MeR),1e-5,'MeR: Moments do not match')
             verifyLessThanOrEqual(testCase,abs(MeH),1e-5,'MeH: Moments do not match')
             
             %unfortunately, it does not seem to be possible to guarantee a
             %higher accuracy. Is this due to numerical error in Hermite
             %interpolation? Or maybe roundoff errors in Gram Schmidt?
             
             % do by hand to verify the results
%              [A,B,C,~,E] = dssdata(sys);
%              
%              V1 = (A-s0(1)*E)\(B*Rt(:,1));
%              V2 = (A-s0(3)*E)\(B*Rt(:,3));
%              V = [real(V1), imag(V1), real(V2), imag(V2)]; V = orth(V);
%              
%              W1 = (A-s0(1)*E).'\(C.'*Lt(:,1));
%              W2 = (A-s0(3)*E).'\(C.'*Lt(:,3));
%              W = [real(W1), imag(W1), real(W2), imag(W2)]; W = orth(W);
%              
%              sysr2 = projectiveMor(sys,V,W);
%              
%              Me2 = moments(sys-sysr2,s0,2);
%              Me = Me2(:,:,1:2:end); Me2(:,:,1:2:end) = [];
%              
%              MeL = mmat(Lt3D,Me);
%              MeR = mmat(Me,Rt3D);
%              MeH = mmat(mmat(Lt3D,Me2),Rt3D);
%              
%              verifyLessThanOrEqual(testCase,abs(MeL),1e-8,'MeL: Moments do not match')
%              verifyLessThanOrEqual(testCase,abs(MeR),1e-8,'MeR: Moments do not match')
%              verifyLessThanOrEqual(testCase,abs(MeH),1e-8,'MeH: Moments do not match')
         end
         function MISO_SIMO (testCase) 
             % test if non-quadratic models are dealt with correctly
             
             sysMIMO = loadSss('CDplayer');
                          
             for idx = 1:2
                 if idx == 1
                     sys = sysMIMO(1:2,1);
                 else
                     sys = sysMIMO(1,1:2);
                 end
             
                 n = 4; 
                 %  Complex shifts and tangential directions
                 Rt = 100*(rand(sys.m,n)+1i*randn(sys.m,n)); Rt(:,[2,4]) = conj(Rt(:,[1,3]));
                 Lt = 100*(rand(sys.p,n)+1i*randn(sys.p,n)); Lt(:,[2,4]) = conj(Lt(:,[1,3]));
                 s0 = 100*(rand(1,n)+1i*randn(1,n)); s0([2,4]) = conj(s0([1,3]));

                 sysr = rk(sys,s0,s0,Rt,Lt);

                 % verify moment matching
                 Me2 = moments(sys-sysr,s0,2);
                 Me = Me2(:,:,1:2:end); Me2(:,:,1:2:end) = [];

                 Lt3D  = permute(Lt.',[3,2,1]);
                 Rt3D  = permute(Rt,[1,3,2]);

                 MeL = mmat(Lt3D,Me);
                 MeR = mmat(Me,Rt3D);
                 MeH = mmat(mmat(Lt3D,Me2),Rt3D);

                 verifyLessThanOrEqual(testCase,abs(MeL),1e-5,'MeL: Moments do not match')
                 verifyLessThanOrEqual(testCase,abs(MeR),1e-5,'MeR: Moments do not match')
                 verifyLessThanOrEqual(testCase,abs(MeH),1e-5,'MeH: Moments do not match')
             end
             end

         function benchmarksSweep2sided (testCase)
             %two-sided reduction for all benchmarks
             for i=1:length(testCase.sysCell)
                 %  test system
                 sys = testCase.sysCell{i};
                 %  get good shifts
                 n = 4; r = ones(sys.m,n); l = ones(sys.p,n);
                 sysrIrka = irka(sys, zeros(1,n),r, l);
                 Opts.rType = 'dir';
                 [r,p] = residue (sysrIrka,Opts);
                 s0 = -(conj(p)); Lt = r{1}; Rt = r{2}.';
                 % make sure real shifts have real directions
                 k = find(imag(s0)==0);
                 if max(imag(Rt(k))) > 1e-10
                     error('Tangential directions corresponding to real shifts are complex!')
                 else
                     Rt(k) = real(Rt(k));
                 end
                 if max(imag(Lt(k))) > 1e-10
                     error('Tangential directions corresponding to real shifts are complex!')
                 else
                     Lt(k) = real(Lt(k));
                 end
                 
                 [sysr, V, W, B_, ~, Rv, C_, ~, Lw] = rk(sys,s0,s0,Rt,Lt);
                 
                 %  Verification with arnoldi
                 [expV,~,~,expW,~] = arnoldi(sys.E,sys.A,sys.B,sys.C,s0,Rt,Lt);
                 
                 % Verification by hand
                 [A,B,C,~,E] = dssdata(sys);

                 V1 = (A-s0(1)*E)\(B*Rt(:,1));
                 V2 = (A-s0(3)*E)\(B*Rt(:,3));
                 Vh = full([real(V1), imag(V1), real(V2), imag(V2)]); Vh = orth(Vh);

                 W1 = (A-s0(1)*E).'\(C.'*Lt(:,1));
                 W2 = (A-s0(3)*E).'\(C.'*Lt(:,3));
                 Wh = full([real(W1), imag(W1), real(W2), imag(W2)]); Wh = orth(Wh);
                 
                 verifyLessThanOrEqual(testCase, subspace(V,Vh),1e-6);
                 verifyLessThanOrEqual(testCase, subspace(W,Wh),1e-6);
                 
                 % The transpose LU problem can be ill conditioned, check the
                 % subspaces instead of the actual matrices!
                 actSolution={sum(svd([V,expV])>1e-12), sum(svd([W,expW])>1e-12)};
                 expSolution={size(V,2), size(W,2)};
                 
                 % Add Sylvester EQ matrices
                 [expRv,expB_] = getSylvester(sys,sysr,V);
                 [expLw,expC_] = getSylvester(sys,sysr,W,'W');
                 
                 actSolution = [actSolution, {Rv, B_, Lw, C_}];
                 expSolution = [expSolution, {expRv,expB_,expLw,expC_}];
                 
                 verification(testCase, actSolution, expSolution, sysr);
                 
                 sysd = sys.'; sysrd = sysr.';
                 res1 = norm(sys.A*V - sys.E*V*(sysr.E\sysr.A) - B_*Rv);
                 res2 = norm(sysd.A*W - sysd.E*W*(sysrd.E\sysrd.A) - C_.'*Lw);
                 
                 verifyEqual(testCase, [res1, res2] , [0, 0], 'AbsTol', 1e-7,...
                     'Sylvester EQ is not satisfied');
                 
                 % verify moment matching only if the "actual" solution
                 % achieves the desired accuracy
                 
                 sysrh = projectiveMor(sys,Vh,Wh);

                 Me2 = moments(sys-sysrh,s0,2);
                 Me = Me2(:,:,1:2:end); Me2(:,:,1:2:end) = [];
                 
                 Lt3D  = permute(Lt.',[3,2,1]);
                 Rt3D  = permute(Rt,[1,3,2]);

                 MeL = mmat(Lt3D,Me);
                 MeR = mmat(Me,Rt3D);
                 MeH = mmat(mmat(Lt3D,Me2),Rt3D);
                 
                 if all(all(abs(MeL)<1e-8)) && all(all(abs(MeR)<1e-8)) && all(all(abs(MeH)<1e-8))
                     Me2 = moments(sys-sysr,s0,2);
                     Me = Me2(:,:,1:2:end); Me2(:,:,1:2:end) = [];

                     MeL = mmat(Lt3D,Me);
                     MeR = mmat(Me,Rt3D);
                     MeH = mmat(mmat(Lt3D,Me2),Rt3D);

                     verifyLessThanOrEqual(testCase,abs(MeL),1e-8,'MeL: Moments do not match')
                     verifyLessThanOrEqual(testCase,abs(MeR),1e-8,'MeR: Moments do not match')
                     verifyLessThanOrEqual(testCase,abs(MeH),1e-8,'MeH: Moments do not match')
                 end
             end
         end          
    end
end

function [] = verification(testCase, actSolution, expSolution, sysr)
       verifyEqual(testCase, actSolution, expSolution,'RelTol',1e-6,'AbsTol',1e-6,...
            'Difference between actual and expected exceeds relative tolerance');
       verifyLessThanOrEqual(testCase, max(imag(sysr.A)), 0, ...
            'Ar is not purely real'); 
       verifyLessThanOrEqual(testCase, max(imag(sysr.E)), 0, ...
            'Er is not purely real'); 
       verifyEqual(testCase, rank(full(sysr.A)), size(sysr.A,1),...
            'Rank(Ar) is not full');
       verifyEqual(testCase, rank(full(sysr.E)), size(sysr.A,1),...
            'Rank(Er) is not full');
       verifyEqual(testCase, nnz(isinf(sysr.A)), 0, ...
            'Ar contains Inf');
       verifyEqual(testCase, nnz(isinf(sysr.E)), 0, ...
            'Er contains Inf');
       verifyEqual(testCase, nnz(isnan(sysr.A)), 0, ...
            'Ar contains Nan');
       verifyEqual(testCase, nnz(isnan(sysr.E)), 0, ...
            'Er contains Nan');
end
function [] = verificationImaginary(testCase, sys, sysr, s0, nMom)
       
       
       % Moment matching
       verifyEqual(testCase, ...
           moments(sys,s0,nMom,Opts),...
           moments(sysr,s0,nMom,Opts),...
           'RelTol',1e-3,'AbsTol',1e-3,...
            sprintf('Moments do not match for %s',sys.Name));
        %reduced order
       verifyEqual(testCase, sysr.n, length(s0)*sysr.m, ... %block krylov
            'Reduced order wrong'); 
        %complex reduced model
       verifyEqual(testCase, all([isreal(sysr.A),isreal(sysr.B),isreal(sysr.C)]), false,...
            'Reduced model is not complex');
end
