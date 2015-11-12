classdef testRk < matlab.unittest.TestCase   
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
%    + test systems: build, beam, random, SpiralInductorPeec 
%      (with E-matrix), LF10 (with E-matrix).
%    + is moment matching achieved
%    + are the matrices of Sylvester equation correct (Bb, Rsylv, Cb, Lsylv)
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
% Last Change:  12 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------
    properties 
        path
        sysCell
    end

    methods(TestClassSetup)
        function getBenchmarks(testCase)
            testCase.path=pwd;
            if exist('benchmarksSysCell.mat','file')
                temp=load('benchmarksSysCell.mat');
                testCase.sysCell=temp.benchmarksSysCell;
            end

            %the directory "benchmark" is in sssMOR
            p = mfilename('fullpath'); k = strfind(p, 'test\'); 
            pathBenchmarks = [p(1:k-1),'benchmarks'];
            cd(pathBenchmarks);
        end
    end
    
    methods(TestClassTeardown)
        function changePath(testCase)
            cd(testCase.path);
        end
    end
    
    methods(Test)
         function testRk1 (testCase) 
              %one-sided reduction isempty(s0_out), s0: real (multiple value)
              load('build.mat'); sys = sss(A,B,C);
              n = 10; s0val = 100; s0 = ones(1,n)*s0val; 
              
              [sysr, V, W, Bb, Rsylv, Cb, Lsylv] = rk(sys, s0);
              actSolution={full(sysr.A), full(sysr.B), full(sysr.C), V, W, ...
                  Rsylv, Bb};
              
              expV = arnoldi(speye(size(A)),A,B, s0);
              [expRsylv,expBb] = getSylvester(sys,sysr,V);
              expSolution={expV'*A*expV, expV'*B, C*expV, expV, expV,...
                   expRsylv, expBb}; 
               
              % check for moment matching as well
              actM = moments(sysr,s0val,n); actSolution = [actSolution,{actM}];
              expM = moments(sys,s0val,n); expSolution = [expSolution,{expM}];
               
              verification(testCase, actSolution, expSolution, sysr);
              verifyEmpty(testCase, Cb, ...
                    'Cb is not empty');
              verifyEmpty(testCase, Lsylv,...
                    'Bt is not empty');
              verifyEqual(testCase, V, W,...
                    'V is not equal W');
         end
         function testRk2 (testCase) 
              %one-sided reduction isempty(s0_in), s0: imag (multiple value)
              load('beam.mat'); sys = sss(A,B,C);
              n = 5; s0val = 100; s0 = [ones(1,n)*s0val*1i,-ones(1,n)*s0val*1i]; 
              
              [sysr, V, W, Bb, Rsylv, Cb, Lsylv] = rk(sys, [], s0);
              actSolution={full(sysr.A), full(sysr.B), full(sysr.C), V, W, Lsylv, Cb};
              
              expW = arnoldi(speye(size(A)),A',C',s0);
              [expLsylv,expCb] = getSylvester(sys,sysr,W,'W');
              expSolution={expW'*A*expW, expW'*B, C*expW, expW, expW, expLsylv,...
                    expCb};
                
              % check for moment matching as well
              actM = moments(sysr,[s0(1),s0(end)], [n,n]); actSolution = [actSolution,{actM}];
              expM = moments(sys,[s0(1),s0(end)], [n,n]); expSolution = [expSolution,{expM}];
                
              verification(testCase, actSolution, expSolution, sysr);
              verifyEmpty(testCase, Bb, ...
                    'Bb is not empty');
              verifyEmpty(testCase, Rsylv,...
                    'Ct is not empty');
              verifyEqual(testCase, V, W,...
                    'V is not equal W');
         end         
         function testRk3 (testCase) 
              %two-sided reduction without E-matrix, s0: zero, imag, real, Inf
              load('random.mat'); sys = sss(A,B,C);
              s0 = [1+4i, 1-4i,4+4i, 4-4i, 0, 0, 5, 300, Inf,  Inf];
              s0moment = s0([1,3,5,7:9]); n = [2, 2, 4, 2, 2, 4];
              
              % careful: Sylvester EQ probably does not hold for shifts at
              %         Infinity!!
              [sysr, V, W] = rk(sys,s0 , s0);
              actSolution={full(sysr.A), full(sysr.B), full(sysr.C), V, W};
              
              [expV,~,expW,~] = arnoldi(speye(size(A)),A,B, C, s0);
              expSolution={expW'*A*expV, expW'*B, C*expV, expV, expW};
              
              % check for moment matching as well
              actM = moments(sysr,s0moment, n); actSolution = [actSolution,{actM}];
              expM = moments(sys,s0moment, n); expSolution = [expSolution,{expM}];

              verification(testCase, actSolution, expSolution, sysr);
              
         end 
         function testRk4 (testCase) 
              %two-sided reduction with E-matrix SISO, s0: zero, imag, real
              load('rail_1357');
              sys = sss(A,B,C,[],E); sys = sys(1,1);
              s0 = [1+4i, 1-4i,4+4i, 4-4i, 0, 0, 5, 300];
              s0moment = s0([1,3,5,7,8]); n = [2, 2, 4, 2, 2];
              
              [sysr, V, W, Bb, Rsylv, Cb, Lsylv] = rk(sys,s0,s0);
              actSolution={full(sysr.A), full(sysr.B), full(sysr.C), full(sysr.E), V, W,...
                  Rsylv, Bb, Lsylv, Cb};
              
              [expV,~,expW,~] = arnoldi(sys.E,sys.A,sys.B,sys.C,s0);
              
              [expRsylv,expBb] = getSylvester(sys,sysr,V);
              [expLsylv,expCb] = getSylvester(sys,sysr,W,'W');
              expSolution={expW'*A*expV, expW'*B, C*expV, expW'*E*expV,expV, expW,...
                  expRsylv,expBb,expLsylv,expCb};
             
              % res = norm(sys.A*V - sys.E*V*(sysr.E\sysr.A)-Bb*Rsylv);
              % res = norm(sys.A*V - sys.E*V*(sysr.E\sysr.A)-expBb*expRsylv);
              
              % check for moment matching as well
              actM = moments(sysr,s0moment, n); actSolution = [actSolution,{actM}];
              expM = moments(sys,s0moment, n); expSolution = [expSolution,{expM}];

              verification(testCase, actSolution, expSolution, sysr);

         end 
%          function testRk5 (testCase) 
              %TEST ERORRED
%               %two-sided reduction with E-matrix, s0: Inf (multiple value)
%               load('SpiralInductorPeec.mat');
% 
%               [expV1,expCt5] = arnoldi(E,A,B,[Inf],@(x,y) (x'*E*y));
%               [expW1,expBt5] = arnoldi(E',A',C',[Inf],@(x,y) (x'*E*y));
% 
%               [sysr, V, W, Bb, Ct, Cb, Bt] = rk(sss(A,B,C,0,E),[Inf] ,[Inf],@(x,y) (x'*E*y));
%               actSolution={W, expW1, full(sysr.A), full(sysr.B), full(sysr.C), full(sysr.E), V, W, Bt, Ct};
%               
% 
%               [expV,expCt5,expW,expBt5] = arnoldi(E,A,B,C,[Inf],@(x,y) (x'*E*y));
%               expSolution={ expW1,expW,expW'*A*expV, expW'*B, C*expV, expW'*E*expV,expV, expW, expBt5,expCt5};
% 
%               verification(testCase, actSolution, expSolution, sysr);
% %               verifyEqual(testCase, {Cb, Bb}, {C - (C*expV)/(expW'*eye(size(A))*expV)*expW',  ...
% %                    B - eye(size(A))*expV*((expW'*eye(size(A))*expV)\(expW'*B))},'AbsTol',0.3,...
% %                    'Cb or Bb failed');
%          end 
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