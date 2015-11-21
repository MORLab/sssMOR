classdef testSs < matlab.unittest.TestCase
    % testSs - testing of ss.m
    %
    % Description:
    %   The function norm.m is tested (3 tests) on:
    %    + Norm of a SISO benchmark system.
    %    + Norm of a SISO random system.
    %    + Norm of a MISO random system.
    %    + Norm of a SIMO random system.
    %    + Norm of MIMO benchmark system.
    %    + Norm of a MIMO random system.
    %    + Verifies for every case the following inputs/outputs: sys=ss(sysSparse)
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
        function testSISObench(testCase)
            load('building.mat');
            sysSparse=sss(A,B,C);
            sys=ss(sysSparse);
        end
        function testSISOrandom(testCase)
            sys=rss(35);
            sysSparse=sss(sys);
             sys=ss(sysSparse);
        end
        function testNormMISOrandom(testCase)
            n=35;
            nInputs=5;
            sys=rss(n);
            sys=ss(sys.A,rand(n,nInputs),sys.C,rand(1,nInputs));
            sysSparse=sss(sys);
             sys=ss(sysSparse);
        end
        function testSIMOrandom(testCase)
            n=35;
            nOutputs=5;
            sys=rss(n);
            sys=ss(sys.A,sys.B,rand(nOutputs,n),rand(nOutputs,1));
            sysSparse=sss(sys);
             sys=ss(sysSparse);
        end
        function testMIMObench(testCase)
            load('cdplayer.mat');
            sysSparse=sss(A,B,C);
             sys=ss(sysSparse);
        end
        function testMIMOrand(testCase)
            n=35;
            nInputs=7;
            nOutputs=5;
            sys=rss(n);
            sys=ss(sys.A,rand(n,nInputs),rand(nOutputs,n),rand(nOutputs,nInputs));
            sysSparse=sss(sys);
            sys=ss(sysSparse);
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution(1:4),  expSolution(1:4),'RelTol',1e-3,...
    'Difference between actual and expected exceeds relative tolerance');
end