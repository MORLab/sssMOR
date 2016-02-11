classdef testTbr < matlab.unittest.TestCase
% testtbr - testing of tbr.m
%
% Description:
%   The function tbr.m is tested (5 tests) on:
%    + q: 5, 10, 15, 20, 25
%    + test systems: building, beam, fom, random, LF10 (with E-matrix).
%    + comparing sysr with directly calculated solution 
%    + W'*V identity matrix
%    + Ar, Er purely real
%    + rank(Ar), rank(Er) full
%    + Neither Inf nor NaN in Ar, Er
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
% Last Change:  11 Feb 2016
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
        function testTbr1(testCase) %q=5
            load('building.mat');      
            q=5;
            
            [sysr, V, W, calchsv] = tbr(sss(A,B,C,0),q);
            load building hsv;
            actSolution={full(sysr.A),full(sysr.B),full(sysr.C),V,W, calchsv};
           
            S=lyapchol(full(A),full(B));
            R=lyapchol(full((A)'),full(C)');
            [Usvd,sigma,Vsvd] = svd(S*R');
            
            expV=S'*Usvd(:,1:q)*(sigma(1:q,1:q))^(-1/2);
            expW=((sigma(1:q,1:q))^(-1/2)*Vsvd(:,1:q)'*R)'; 
            
            expSolution={expW'*full(A)*expV, expW'*full(B),full(C)*expV,  expV, expW, hsv};
            verification(testCase, actSolution, expSolution, sysr);
            verifyLessThan(testCase, norm((W'*V)-eye(size(V,2))),1.01,...
               'W*V not identity matrix');
        end
        
        function testTbr2(testCase) %q=15
            load('beam.mat');
            q=15;
            
            [sysr, V, W] = tbr(sss(A,B,C,0),q);
            actSolution={full(sysr.A),full(sysr.B),full(sysr.C),V,W};

            S=lyapchol(full(A),full(B));
            R=lyapchol(full((A)'),full(C)');
            [Usvd,sigma,Vsvd] = svd(S*R');
            
            expV=S'*Usvd(:,1:q)*(sigma(1:q,1:q))^(-1/2);
            expW=((sigma(1:q,1:q))^(-1/2)*Vsvd(:,1:q)'*R)'; 
            
            expSolution={expW'*full(A)*expV, expW'*full(B),full(C)*expV,  expV, expW};
            verification(testCase, actSolution, expSolution, sysr);
            verifyLessThan(testCase, norm((W'*V)-eye(size(V,2))),1.01,...
               'W*V not identity matrix');
        end        
        function testTbr3(testCase) %q=25
            load('fom.mat');
            q=25;
            Opts.adi=0;
            
            [sysr, V, W] = tbr(sss(A,B,C,0),q, Opts);
            actSolution={full(sysr.A),full(sysr.B),full(sysr.C),V,W};
            
            S=lyapchol(full(A),full(B));
            R=lyapchol(full((A)'),full(C)');
            [Usvd,sigma,Vsvd] = svd(S*R');
            
            expV=S'*Usvd(:,1:q)*(sigma(1:q,1:q))^(-1/2);
            expW=((sigma(1:q,1:q))^(-1/2)*Vsvd(:,1:q)'*R)'; 
            
            expSolution={expW'*full(A)*expV, expW'*full(B),full(C)*expV,  expV, expW};
            verification(testCase, actSolution, expSolution, sysr);
            verifyLessThan(testCase, norm((W'*V)-eye(size(V,2))),1.01,...
                'W*V not identity matrix');
        end
        function testTbr4(testCase) %q=20
            load('random.mat');
            q=20;
            
            [sysr, V, W] = tbr(sss(A,B,C,0),q);
            actSolution={full(sysr.A),full(sysr.B),full(sysr.C),V,W};
            
            S=lyapchol(full(A),full(B));
            R=lyapchol(full((A)'),full(C)');
            [Usvd,sigma,Vsvd] = svd(S*R');
            
            expV=S'*Usvd(:,1:q)*(sigma(1:q,1:q))^(-1/2);
            expW=((sigma(1:q,1:q))^(-1/2)*Vsvd(:,1:q)'*R)'; 

            expSolution={expW'*full(A)*expV, expW'*full(B),full(C)*expV,  expV, expW};
            verification(testCase, actSolution, expSolution, sysr);
            verifyLessThan(testCase, norm((W'*V)-eye(size(V,2))),1.01,...
               'W*V not identity matrix');
        end
        function testTbr5(testCase) %q=10 (with E-matrix)
            load('LF10.mat');
            q=10;
            
            E=blkdiag(speye(size(M)), M);
            A=[zeros(size(M)),speye(size(M)); -K, -D];
            B=[zeros(size(M,1),1); B];
            C=[C, zeros(1,size(M,1))];
            
            [sysr, V, W] = tbr(sss(E\A,E\B,C,0),q);
            actSolution={full(sysr.A),full(sysr.B),full(sysr.C),V,(W'/E)'};
            
            S=lyapchol(full(E\A),full(E\B));
            R=lyapchol(full((E\A)'),full(C)');
            [Usvd,sigma,Vsvd] = svd(S*R');
            
            expV=S'*Usvd(:,1:q)*(sigma(1:q,1:q))^(-1/2);
            expW=((sigma(1:q,1:q))^(-1/2)*Vsvd(:,1:q)'*R/E)'; %warum W/E?
            
            expSolution={expW'*full(A)*expV, expW'*full(B),full(C)*expV, expV, expW};
            
            verification(testCase, actSolution, expSolution, sysr);
            verifyEqual(testCase, full(sysr.E),  expW'*E*expV,'RelTol',0.2,'AbsTol',0.0000000001,...
                    'sysr.E');
        end
        function testTbr6(testCase) %adi
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae && sys.n>100
                    q=50;
                    opts.adi='adi';
                    [~,~,~,actHsv]=tbr(sys,q,opts);
                    opts.adi=0;
                    [~,~,~,expHsv]=tbr(sys,q,opts);
                    
                    actSolution={actHsv(1:5)};
                    expSolution={expHsv(1:5)};

                    verifyEqual(testCase, actSolution, expSolution,'RelTol',0.3,...
                        'Difference between actual and expected exceeds relative tolerance');
                end
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution, sysr)
       verifyEqual(testCase, actSolution, expSolution,'RelTol',0.3,...
            'Difference between actual and expected exceeds relative tolerance');
       verifyLessThanOrEqual(testCase, max(imag(sysr.A)), 0, ...
            'Ar is not purely real'); 
       verifyLessThanOrEqual(testCase, max(imag(sysr.E)), 0, ...
            'Er is not purely real'); 
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