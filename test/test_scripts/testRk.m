classdef testRk < matlab.unittest.TestCase   
% testRk - testing of rk.m
%
% Description:
%   The function rk.m is tested (5 tests) on:
%    + comparing sysr to a directly calculated solution based on arnoldi.m
%    + TO DO: Ct, Bt, Cb, Bb (?)
%    + one-sided (V=W), double-sided 
%    + Ar, Er purely real
%    + rank(Ar), rank(Er) full
%    + Neither Inf nor NaN in Ar, Er
%    + s0: purely real, purely imaginary, zero, Inf 
%    + test systems: build, beam, random, SpiralInductorPeec 
%      (with E-matrix), LF10 (with E-matrix).
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
% Last Change:  05 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------
    properties
    end
 
 
    methods(Test)
         function testRk1 (testCase) 
              %one-sided reduction isempty(s0_out), s0: real (multiple value)
              load('build.mat');
              [sysr, V, W, Bb, Ct, Cb, Bt] = rk(sss(A,B,C), ones(1,10)*100,[]);
              actSolution={full(sysr.A), full(sysr.B), full(sysr.C), V, W, Ct, Bb};
              
              [expV,expCt] = arnoldi(speye(size(A)),A,B,ones(1,10)*100);
              expSolution={expV'*A*expV, expV'*B, C*expV, expV, expV, expCt,...
                   B - speye(size(A))*expV*(speye(size(sysr.A))\sysr.B)}; 
              verification(testCase, actSolution, expSolution, sysr);
              verifyEmpty(testCase, Cb, ...
                    'Cb is not empty');
              verifyEmpty(testCase, Bt,...
                    'Bt is not empty');
              verifyEqual(testCase, V, W,...
                    'V is not equal W');
         end
         function testRk2 (testCase) 
              %one-sided reduction isempty(s0_in), s0: imag (multiple value)
              load('beam.mat');
              [sysr, V, W, Bb, Ct, Cb, Bt] = rk(sss(A,B,C), [], [ones(1,5)*100*1i,-ones(1,5)*100*1i]);
              actSolution={full(sysr.A), full(sysr.B), full(sysr.C), V, W, Bt, Cb};
              
              [expW,expBt2] = arnoldi(speye(size(A)),A',C',[ones(1,5)*100*1i,-ones(1,5)*100*1i]);
               expSolution={expW'*A*expW, expW'*B, C*expW, expW, expW, expBt2',...
                    C - (C*expW)/(expW'*eye(size(A))*expW)*expW'};
              verification(testCase, actSolution, expSolution, sysr);
              verifyEmpty(testCase, Bb, ...
                    'Bb is not empty');
              verifyEmpty(testCase, Ct,...
                    'Ct is not empty');
              verifyEqual(testCase, V, W,...
                    'V is not equal W');
         end
         
         function testRk3 (testCase) 
              %two-sided reduction without E-matrix, s0: zero, imag, real, Inf
              load('random.mat');
              [sysr, V, W, Bb, Ct, Cb, Bt] = rk(sss(A,B,C),[1+4i, 1-4i,4+4i, 4-4i, 0, 0, 5, 300, Inf,  Inf] ,...
                   [1+4i, 1-4i,4+4i, 4-4i, 0, 0, 5, 300, Inf,  Inf]);
              actSolution={full(sysr.A), full(sysr.B), full(sysr.C), V, W, Bt, Ct};
              
%               [expV,~,~,expCt5] = arnoldi(speye(size(A)),A,B,[1+4i, 1-4i,4+4i, 4-4i, 0, 0, 5, 300, Inf,  Inf]);
              [expV,expCt5,expW,expBt5] = arnoldi(speye(size(A)),A,B, C,[1+4i, 1-4i,4+4i, 4-4i, 0, 0, 5, 300, Inf,  Inf]);
              expSolution={expW'*A*expV, expW'*B, C*expV, expV, expW, expBt5,expCt5};

              verification(testCase, actSolution, expSolution, sysr);
              verifyEqual(testCase, {Cb, Bb}, {C - (C*expV)/(expW'*eye(size(A))*expV)*expW',  ...
                   B - eye(size(A))*expV*((expW'*eye(size(A))*expV)\(expW'*B))},'AbsTol',0.000001,'RelTol',0.3,...
                   'Cb or Bb failed');
         end 
         function testRk4 (testCase) 
              %two-sided reduction with E-matrix, s0: zero, imag, real
              load('LF10.mat');
              E=blkdiag(speye(size(M)), M);
              A=[zeros(size(M)),speye(size(M)); -K, -D];
              B=[zeros(size(M,1),1); B];
              C=[C, zeros(1,size(M,1))];
              
              [sysr, V, W, Bb, Ct, Cb, Bt] = rk(sss(A,B,C,0,E),[1+4i, 1-4i,4+4i, 4-4i, 0, 0, 5, 300] ,...
                   [1+4i, 1-4i,4+4i, 4-4i, 0, 0, 5, 300]);
              actSolution={full(sysr.A), full(sysr.B), full(sysr.C), full(sysr.E), V, W, Bt, Ct};
              
              [expV,expCt5,expW,expBt5] = arnoldi(E,A,B,C,[1+4i, 1-4i,4+4i, 4-4i, 0, 0, 5, 300]);
              expSolution={expW'*A*expV, expW'*B, C*expV, expW'*E*expV,expV, expW, expBt5,expCt5};

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
       verifyEqual(testCase, actSolution, expSolution,'RelTol',0.2,'AbsTol',0.00000001,...
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