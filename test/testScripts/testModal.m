classdef testModal < matlab.unittest.TestCase
% testModal - testing of modalMor.m
%
% Description:
%   The function modalMor.m is tested (3 tests) on:
%    + comparing the eigenvalues of the reduced system to the solution of
%      modreal (only 'SM' possible).
%    + test systems: diagonal, building, LF10 (with E-matrix)
%
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
% Last Change:  07 Sep 2015
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
            if isempty(testCase.sysCell)
                error('No benchmarks loaded.');
            end

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
        function testModal1(testCase) 
            %Diagonal matrix
            A=zeros(11);
            A=diag(-10:-10:-100);
            A(11,11)=-30;
            B=(1:11)';
            C=B';

            Opts.type='SM';
            [sysr] = modalMor(sss(A,B,C,0), 6, Opts);
            actSolution={full(sort(eig(sysr)))};
            
            [expsysr,~]=modreal(ss(full(A),full(B),full(C),0),6);
            expSolution={full(sort(eig(expsysr)))};
                     
            verification(testCase, actSolution, expSolution, sysr);
        end
        function testModal2(testCase) 
            %without E-matrix
            load('building.mat');

            Opts.type='SM';
            [sysr] = modalMor(sss(A,B,C,0), 6, Opts);
            actEig=sort(eig(sysr));
            actSolution={full(real(actEig)), full(abs(imag(actEig)))};
            
            [expsysr,~]=modreal(ss(full(A),full(B),full(C),0),6);
            expEig=sort(eig(expsysr));
            expSolution={full(real(expEig)), full(abs(imag(expEig)))};
                 
            verification(testCase, actSolution, expSolution, sysr);
        end
        function testModal3(testCase) 
            %with E-matrix
            load('LF10.mat');
            E=blkdiag(speye(size(M)), M);
            A=[zeros(size(M)),speye(size(M)); -K, -D];
            B=[zeros(size(M,1),1); B];
            C=[C, zeros(1,size(M,1))];
            
            Opts.type='SM';
            [sysr] = modalMor(sss(A,B,C,0,E), 9, Opts);
            actSolution=full(sort(eig(sysr)));
%             actSolution={real(actEig), abs(imag(actEig))};
            
            [expsysr,~]=modreal(ss(full(E\A),full(E\B),full(C),0),9);
            expSolution=full(sort(eig(expsysr)));
%             expSolution={real(expEig), abs(imag(expEig))};
                 
            verification(testCase, actSolution, expSolution, sysr);
        end
    end
end

function [] = verification(testCase, actSolution, expSolution, sysr)
       verifyEqual(testCase, actSolution, expSolution,'RelTol',0.2,'AbsTol',0.00000001,...
            'Difference between actual and expected exceeds relative tolerance');
       verifyEqual(testCase, rank(full(sysr.A)), length(sysr.B),...
            'Rank(Ar) is not full');
       verifyEqual(testCase, rank(full(sysr.E)), length(sysr.B),...
            'Rank(Er) is not full');
       verifyEqual(testCase, nnz(isinf(sysr.A)), 0, ...
            'Ar contains Inf');
       verifyEqual(testCase, nnz(isinf(sysr.E)), 0, ...
            'Er contains Inf');
       verifyEqual(testCase, nnz(isnan(sysr.A)), 0, ...
            'Ar contains Nan');
       verifyEqual(testCase, nnz(isnan(sysr.E)), 0, ...
            'Er contains Nan');
       verifyLessThan(testCase, max(real(eig(sysr))), 0, ...
            'Ar is not stable');
end