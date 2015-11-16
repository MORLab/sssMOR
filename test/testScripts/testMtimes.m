classdef testMtimes < matlab.unittest.TestCase
    % testMtimes - testing of mtimes.m
%
% Description:
%   The function rk.m is tested (3 tests) on:
%    + combination of two benchmark-systems.
%    + combination of two random-systems that are equal.
%    + combination of two random-systems that are different.
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
%               Jorge Luiz Moreira Silva
% Last Change:  26 Out 2015
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
                temp=load('benchmarksSysCell.mat');
                testCase.sysCell=temp.benchmarksSysCell;
                testCase.deleteBenchmarks=0;
            else
                testCase.testPath=loadBenchmarks;
                testCase.deleteBenchmarks=1;
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
        function testMTimes1(testCase)
            load('build.mat');
            sysSparse=sss(A,B,C);
            sys=ss(full(A),full(B),full(C),zeros(1,1));
            
            resultSparse = mtimes(sysSparse, sysSparse);
            result=mtimes(sys,sys);
            verification(testCase, resultSparse, result);
        end
        function testMTimes2(testCase)
            sys=rss(50);
            sysSparse=sss(sys);
            
            resultSparse = mtimes(sysSparse, sysSparse);
            result=mtimes(sys,sys);
            verification(testCase, resultSparse, result);
        end
        function testMTimes3(testCase)
            sys1=rss(50);
            sys2=rss(35);
            sysSparse1=sss(sys1);
            sysSparse2=sss(sys2);
            
            resultSparse = mtimes(sysSparse1, sysSparse2);
            result=mtimes(sys1,sys2);
            verification(testCase, resultSparse, result);
        end
    end
end

function [] = verification(testCase, actSolution, expSolution, m)
verifyEqual(testCase, full(actSolution.A),  full(expSolution.A),'RelTol',0.1,...
    'Difference between actual and expected exceeds relative tolerance');
verifyEqual(testCase,  full(actSolution.B),  full(expSolution.B),'RelTol',0.1,...
    'Difference between actual and expected exceeds relative tolerance');
verifyEqual(testCase,  full(actSolution.C),  full(expSolution.C),'RelTol',0.1,...
    'Difference between actual and expected exceeds relative tolerance');
verifyEqual(testCase,  full(actSolution.D),  full(expSolution.D),'RelTol',0.1,...
    'Difference between actual and expected exceeds relative tolerance');
end