classdef testIrka < matlab.unittest.TestCase
% testIrka - testing of irka.m
%
% Description:
%   The function irka.m is tested (5 tests) on:
%    + comparing sysr to a directly calculated solution based on RK.m and
%      the final s0
%    + k>kmax, error>epsilon
%    + Ar, Er purely real
%    + rank(Ar), rank(Er) full
%    + Neither Inf nor NaN in Ar, Er
%    + s0: purely real, purely imaginary, zero, Inf (Markov-parameter)
%    + test systems: build, beam, random, SpiralInductorPeec 
%      (includes E-matrix), LF10 (includes E-matrix).
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
% Last Change:  05 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------ 

    properties 
        pwdPath
        sysCell
        deleteBenchmarks
        testPath
    end

    methods(TestClassSetup)
        function getBenchmarks(testCase)
            testCase.pwdPath=pwd;
            if exist('benchmarksSysCell.mat','file')
                testCase.deleteBenchmarks=0;
            else
                testCase.testPath=loadBenchmarks;
                testCase.deleteBenchmarks=1;
            end
            
            temp=load('benchmarksSysCell.mat');
            testCase.sysCell=temp.benchmarksSysCell;

            %the directory "benchmark" is in sssMOR
            p = mfilename('fullpath'); k = strfind(p, 'test\'); 
            pathBenchmarks = [p(1:k-1),'benchmarks'];
            cd(pathBenchmarks);
        end
    end
    
    methods(TestClassTeardown)
        function changePath(testCase)
            if testCase.deleteBenchmarks
                cd(testCase.testPath);
                delete('benchmarksSysCell.mat');
            end
            cd(testCase.pwdPath);
        end
    end
 
    methods(Test)
        function testIrka1(testCase) 
             %s0: zero, real, imag, Inf
            load('build.mat');
            Opts=struct('epsilon',0.05,'maxiter',300,'type','stab','stopCrit','s0');

            [sysr, V, W, s0, s0_traj] = irka(sss(A,B,C), ...
                 [0,0,50,100,Inf, Inf,1+5i,1-5i,14-0.2i,14+0.2i],Opts);
            actSolution={full(sysr.A),full(sysr.B),full(sysr.C),V,W};
            
            expV=arnoldi(eye(size(A)),A,B,s0);
            expW=arnoldi(eye(size(A)),A',C',s0);
            expSolution={expW'*A*expV, expW'*B, C*expV, expV, expW};
            
            verification(testCase, actSolution, expSolution, sysr, s0, s0_traj);
            verifyLessThan(testCase, max(real(eig(sysr))), 0, ...
                    'Ar is not stable');
        end
        function testIrka2(testCase) 
             %s0: Inf (multiple value), imag (multiple value), real (multiple value)
            load('random.mat');
            opts=struct('epsilon',0.05,'maxiter',300,'type','stab','stopCrit','s0');

            [sysr, V, W, s0, s0_traj] = irka(sss(A,B,C), [Inf, Inf, 3+5i, 3-5i, 3+5i, 3-5i, 60, 60, 60],opts);
            actSolution={full(sysr.A),full(sysr.B),full(sysr.C),V,W};
            
            expV=arnoldi(eye(size(A)),A,B,s0);
            expW=arnoldi(eye(size(A)),A',C',s0);
            expSolution={expW'*A*expV, expW'*B, C*expV, expV, expW};
            
            verification(testCase, actSolution, expSolution, sysr, s0, s0_traj);
            verifyLessThan(testCase, max(real(eig(sysr))), 0, ...
                    'Ar is not stable');
        end
        function testIrka3(testCase) 
             %s0: real, imag, zero, Inf
            load('fom.mat');
            opts=struct('epsilon',0.05,'maxiter',300,'type','orig','stopCrit','s0');
               
            [sysr, V, W, s0, s0_traj] = irka(sss(A,B,C), ...
                 [0,0,3,300,Inf, Inf, 10+50i,10-50i], opts);
             
            expV=arnoldi(eye(size(A)),A,B,s0);
            expW=arnoldi(eye(size(A)),A',C',s0);
            
            actSolution={rank([V,expV]),rank([W,expW])};
            expSolution={size(V,2), size(W,2)};
            
            verification(testCase, actSolution, expSolution, sysr, s0, s0_traj);
        end
        function testIrka4(testCase) 
             %s0: real, imag, Inf, zero (with E-matrix)
            load('SpiralInductorPeec.mat');
            opts=struct('epsilon',0.05,'maxiter',300,'type','orig','stopCrit','s0');

            [sysr, V, W, s0, s0_traj] = irka(sss(A, B, C,0,E), ...
                 [0,300000,25603+546i,25603-546i,Inf], opts);

            expV=arnoldi(E,A,B,s0);
            expW=arnoldi(E',A',C',s0);
            svdV = svd([V,expV],0); svdW = svd([W,expW]); tol = 1e-6;
            rV = sum(svdV>tol); rW = sum(svdW>tol);
            actSolution={rV,rW};
            expSolution={size(V,2), size(W,2)};
            
            verification(testCase, actSolution, expSolution, sysr, s0, s0_traj);
        end
        function testIrka5(testCase) 
             %s0: Inf (multiple value), imag (multiple value), real (multiple value)
            load('LF10.mat');
            opts=struct('epsilon',0.05,'maxiter',300,'type','stab','stopCrit','s0');
            
            E=blkdiag(speye(size(M)), M);
            A=[zeros(size(M)),speye(size(M)); -K, -D];
            B=[zeros(size(M,1),1); B];
            C=[C, zeros(1,size(M,1))];

            [sysr, V, W, s0, s0_traj] = irka(sss(A, B, C,0,E), ...
                 [Inf, Inf, 3+5i, 3-5i, 3+5i, 3-5i, 60, 60, 60], opts);
            
            expV=arnoldi(E,A,B,s0);
            expW=arnoldi(E',A',C',s0);
            
            svdV = svd([V,expV],0); svdW = svd([W,expW]); tol = 1e-6;
            rV = sum(svdV>tol); rW = sum(svdW>tol);
            actSolution={rV,rW};
            expSolution={size(V,2), size(W,2)};
            
            verification(testCase, actSolution, expSolution, sysr, s0, s0_traj);
        end
        function testIrka6(testCase)
            % all benchmarks + check hermite interpolation
            for i=1:length(testCase.sysCell)
                %  test system
                sys=testCase.sysCell{i};
                
                Opts=struct('epsilon',0.05,'maxiter',300,'type','stab','stopCrit','s0');
                [sysr, V, W, s0, s0_traj, Rt, Lt] = irka(sys,zeros(1,10),ones(sys.m,10), ones(sys.p,10),Opts);
                
                expV=arnoldi(sys.E,sys.A,sys.B,s0, Rt); 
                expW=arnoldi(sys.E',sys.A',sys.C',s0, Lt);
                            
                svdV = svd([V,expV],0); svdW = svd([W,expW]); tol = 1e-6;
                rV = sum(svdV>tol); rW = sum(svdW>tol);
                actSolution={rV,rW};
                expSolution={size(V,2), size(W,2)};
              
              % verify moment matching
              actM = moments(sysr,s0, 2); expM = moments(sys,s0, 2);
              actMt = {}; expMt = {};
              for iS = 1:length(s0)
                  jM = iS*2-1;
                  actMt = [actMt, {actM(:,:,jM)*Rt(:,iS), Lt(:,iS).'*actM(:,:,jM),...
                           Lt(:,iS).'*actM(:,:,jM)*Rt(:,iS), Lt(:,iS).'*actM(:,:,jM+1)*Rt(:,iS)}];
                  expMt = [expMt, {expM(:,:,jM)*Rt(:,iS), Lt(:,iS).'*expM(:,:,jM),...
                           Lt(:,iS).'*expM(:,:,jM)*Rt(:,iS), Lt(:,iS).'*expM(:,:,jM+1)*Rt(:,iS)}];
              end
              actSolution = [actSolution, actMt];
              expSolution = [expSolution, expMt];
              verification(testCase, actSolution, expSolution, sysr, s0, s0_traj);
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution, sysr, s0, s0_traj)
       verifyEqual(testCase, actSolution, expSolution,'RelTol',0.3,'AbsTol',0.00000001,...
            'Difference between actual and expected exceeds relative tolerance');
       verifyLessThanOrEqual(testCase, max(imag(sysr.A)), 0, ...
            'Ar is not purely real'); 
       verifyLessThanOrEqual(testCase, max(imag(sysr.E)), 0, ...
            'Er is not purely real'); 
       verifyEqual(testCase, rank(full(sysr.A)), length(s0),...
            'Rank(Ar) is not full');
       verifyEqual(testCase, rank(full(sysr.E)), length(s0),...
            'Rank(Er) is not full');
       verifyEqual(testCase, nnz(isinf(sysr.A)), 0, ...
            'Ar contains Inf');
       verifyEqual(testCase, nnz(isinf(sysr.E)), 0, ...
            'Er contains Inf');
       verifyEqual(testCase, nnz(isnan(sysr.A)), 0, ...
            'Ar contains Nan');
       verifyEqual(testCase, nnz(isnan(sysr.E)), 0, ...
            'Er contains Nan');
       verifyLessThanOrEqual(testCase, size(s0_traj,1),301, ...
            'k > kmax');
       verifyLessThanOrEqual(testCase, norm((s0_traj(end,:)-s0)./s0_traj(end,:), 1)/length(s0), 0.05,...
            'error > epsilon - irka failed to converge');
end