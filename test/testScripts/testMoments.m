classdef testMoments < matlab.unittest.TestCase
% testMoments - testing of moments.m using all benchmarks in local folder
%
% Description:
%   The function moments.m is tested (4 tests) on:
%    + moments: mi=c*(A-s0*I)^-(i+1)*b i=0,1,...
%    + markov: mi=c*A^i*b, i=0,1,..
%    + test systems: build, beam, random, LF10 (with E-matrix)
%    + Neither Inf nor NaN in m
%    + s0: real, imag, Inf
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto, Lisa Jeschek
% Last Change:  31 Oct 2015
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
        function testMoments1(testCase) 
            for i=1:length(testCase.sysCell)  
                if testCase.sysCell{i}.isDescriptor==0
                    sys=testCase.sysCell{i};
                    s0=5;

                    m = moments(sys, s0, 4);
                    actSolution={m};
                    
                    expm = expmoments(sys,s0,4);
                    expSolution={expm};
                    
                    verification(testCase, actSolution, expSolution, m);
                end
            end
        end
        
        function testMoments2(testCase) 
            %s0: imag
            for i=1:length(testCase.sysCell)  
                if testCase.sysCell{i}.isDescriptor==0
                    sys=testCase.sysCell{i};  
                    s0=5+7i;

                    m = moments(sys, s0, 4);
                    actSolution={m};

                    expm = expmoments(sys,s0,4);
                    expSolution={expm};
                    verification(testCase, actSolution, expSolution, m);
                end
            end
        end
        
        function testMoments3(testCase) 
            %s0: Inf
            for i=1:length(testCase.sysCell)  
                if testCase.sysCell{i}.isDescriptor==0
                    sys=testCase.sysCell{i};    
                    s0=Inf;

                    m = moments(sys, s0, 4);
                    actSolution={m};
                    
                    
                    expm = expmoments(sys,s0,4);
                    expSolution={expm};
                    verification(testCase, actSolution, expSolution, m);
                end
            end
        end
        
         function testMoments4(testCase) 
            %with E-matrix, s0 real
            for i=1:length(testCase.sysCell)  
                if testCase.sysCell{i}.isDescriptor==1
                    sys=testCase.sysCell{i};    
                    s0=50;

                    m = moments(sys, s0, 4);
                    actSolution={m};

                    expm = expmoments(sys, s0, 4);
                    expSolution={expm};
                    verification(testCase, actSolution, expSolution, m);
                end  
            end
         end
    end 
end

function [] = verification(testCase, actSolution, expSolution, m)
       verifyEqual(testCase, actSolution, expSolution,'RelTol',0.1,...
            'Difference between actual and expected exceeds relative tolerance');
       verifyEqual(testCase, nnz(isinf(m)), 0, ...
            'm contains Inf');
       verifyEqual(testCase, nnz(isnan(m)), 0, ...
            'm contains Nan');
end

function expm = expmoments(sys,s0,n)
        expm=zeros(size(sys.C,1),size(sys.B,2),n);        
        if isinf(s0) && ~sys.isDescriptor
            %mi=c*A^i*b
            expm(:,:,1)=0; 
            expm(:,:,2)=sys.C*sys.B;
            expm(:,:,3)=sys.C*sys.A*sys.B;
            expm(:,:,4)=sys.C*sys.A^2*sys.B;
        else
            [L,U] = lu(sys.A-s0*sys.E); temp = (U\(L\sys.B));
            for iMoment = 1:n
                expm(:,:,iMoment)=sys.C*temp;
                temp = (U\(L\(sys.E*temp)));
            end
        end
end