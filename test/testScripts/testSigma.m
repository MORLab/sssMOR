classdef testSigma < matlab.unittest.TestCase
    % testSigma - testing of signa.m
    %
    % Description:
    %   The function sigma.m is tested (3 tests) on:
    %    + Norm of a SISO benchmark system.
    %    + Norm of a SISO random system.
    %    + Norm of a MISO random system.
    %    + Norm of a SIMO random system.
    %    + Norm of MIMO benchmark system.
    %    + Norm of a MIMO random system.
    %    + Verifies for every case just if it runs for the syntax
    %    [mag,omega]=sigma(sys) and sigma(sys)
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
         sysCell;
         testPath;
    end
 
    methods(TestMethodSetup)
        function getBenchmarks(testCase)
            testCase.testPath=pwd;
            if exist('benchmarksSysCell.mat','file')
                load('benchmarksSysCell.mat');
                testCase.sysCell=benchmarksSysCell;
            end
            
            %the directory "benchmark" is in sssMOR
            p = mfilename('fullpath'); k = strfind(p, 'test\'); 
            pathBenchmarks = [p(1:k-1),'benchmarks'];
            cd(pathBenchmarks);
        end
    end
    
    methods(TestMethodTeardown)
        function changePath(testCase)
            cd(testCase.testPath);
        end
    end
    
    methods(Test)
        function testSISObench(testCase)
            load('build.mat');
            sysSparse=sss(A,B,C);
            sys=ss(A,B,C,zeros(1,1));
            sigma(sysSparse);
            [mag,omega]=sigma(sysSparse);
            close all;
        end
        function testSISOrandom(testCase)
            sys=rss(35);
            sysSparse=sss(sys);
            sigma(sysSparse);
            [mag,omega]=sigma(sysSparse);
            close all;
        end
        function testMISO(testCase)
            n=35;
            nInputs=5;
            sys=rss(n);
            sys=ss(sys.A,rand(n,nInputs),sys.C,rand(1,nInputs));
            sysSparse=sss(sys);
            sigma(sysSparse);
            [mag,omega]=sigma(sysSparse);
            close all;
        end
        function testSIMO(testCase)
            n=35;
            nOutputs=5;
            sys=rss(n);
            sys=ss(sys.A,sys.B,rand(nOutputs,n),rand(nOutputs,1));
            sysSparse=sss(sys);
            sigma(sysSparse);
            [mag,omega]=sigma(sysSparse);
            close all;
        end
        function testMIMObench(testCase)
            load('cdplayer.mat');
            sysSparse=sss(A,B,C);
            sys=ss(full(A),full(B),full(C),zeros(2,2));
            sigma(sysSparse);
            [mag,omega]=sigma(sysSparse);
            close all;
        end
        function testMIMOrandom(testCase)
            n=35;
            nInputs=7;
            nOutputs=5;
            sys=rss(n);
            sys=ss(sys.A,rand(n,nInputs),rand(nOutputs,n),rand(nOutputs,nInputs));
            sysSparse=sss(sys);
            sigma(sysSparse);
            [mag,omega]=sigma(sysSparse);
            close all;
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution(1:4),  expSolution(1:4),'RelTol',1e-3,...
    'Difference between actual and expected exceeds relative tolerance');
end