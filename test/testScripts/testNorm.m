classdef testNorm < matlab.unittest.TestCase
    % testNorm - testing of norm.m
    %
    % Description:
    %   The function norm.m is tested (3 tests) on:
    %    + Norm of a SISO benchmark system.
    %    + Norm of a SISO random system.
    %    + Norm of a MISO random system.
    %    + Norm of a SIMO random system.
    %    + Norm of MIMO benchmark system.
    %    + Norm of a MIMO random system.
    %    + Verifies for every case the following inputs/outputs (norm(sys),
    %       norm(sys,inf), [n,fpeak]=norm(sys,inf), norm(sys,2)
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
            sys=ss(A,B,C,zeros(1,1));
            actNorm1=norm(sysSparse);
            [actNorm2,actFreq]=norm(sysSparse,inf);
            sysSparse.D=zeros(size(sysSparse.D));
            actNorm3=norm(sysSparse,2);
            actNorm = [actNorm1,actNorm2,actFreq,actNorm3];
            expNorm1=norm(sys);
            [expNorm2,expFreq]=norm(sys,inf);
            sys.D=zeros(size(sys.D));
            expNorm3=norm(sys,2);
            if norm(freqresp(sys,actFreq),2)>expNorm2
                expNorm2=actNorm2;
                expFreq=actFreq;
            end
            if (actFreq==0 &&expFreq==0)
                expNorm2=actNorm2;
            end
            expNorm = [expNorm1,expNorm2,expFreq,expNorm3];
            verification(testCase, actNorm, expNorm);
        end
        function testSISOrandom(testCase)
            sys=rss(35);
            sysSparse=sss(sys);
            actNorm1=norm(sysSparse);
            [actNorm2,actFreq]=norm(sysSparse,inf);
            sysSparse.D=zeros(size(sysSparse.D));
            actNorm3=norm(sysSparse,2);
            actNorm = [actNorm1,actNorm2,actFreq,actNorm3];
            expNorm1=norm(sys);
            [expNorm2,expFreq]=norm(sys,inf);
            sys.D=zeros(size(sys.D));
            expNorm3=norm(sys,2);
            if norm(freqresp(sys,actFreq),2)>expNorm2
                expNorm2=actNorm2;
                expFreq=actFreq;
            end
            if (actFreq==0 &&expFreq==0)
                expNorm2=actNorm2;
            end
            expNorm = [expNorm1,expNorm2,expFreq,expNorm3];
            verification(testCase, actNorm, expNorm);
        end
        function testMISOrandom(testCase)
            n=35;
            nInputs=5;
            sys=rss(n);
            sys=ss(sys.A,rand(n,nInputs),sys.C,rand(1,nInputs));
            sysSparse=sss(sys);
            actNorm1=norm(sysSparse);
            [actNorm2,actFreq]=norm(sysSparse,inf);
            sysSparse.D=zeros(size(sysSparse.D));
            actNorm3=norm(sysSparse,2);
            actNorm = [actNorm1,actNorm2,actFreq,actNorm3];
            expNorm1=norm(sys);
            [expNorm2,expFreq]=norm(sys,inf);
            sys.D=zeros(size(sys.D));
            expNorm3=norm(sys,2);
            if norm(freqresp(sys,actFreq),2)>expNorm2
                expNorm2=actNorm2;
                expFreq=actFreq;
            end
            if (actFreq==0 &&expFreq==0)
                expNorm2=actNorm2;
            end
            expNorm = [expNorm1,expNorm2,expFreq,expNorm3];
            verification(testCase, actNorm, expNorm);
        end
        function testSIMOrandom(testCase)
            n=35;
            nOutputs=5;
            sys=rss(n);
            sys=ss(sys.A,sys.B,rand(nOutputs,n),rand(nOutputs,1));
            sysSparse=sss(sys);
            actNorm1=norm(sysSparse);
            [actNorm2,actFreq]=norm(sysSparse,inf);
            sysSparse.D=zeros(size(sysSparse.D));
            actNorm3=norm(sysSparse,2);
            actNorm = [actNorm1,actNorm2,actFreq,actNorm3];
            expNorm1=norm(sys);
            [expNorm2,expFreq]=norm(sys,inf);
            sys.D=zeros(size(sys.D));
            expNorm3=norm(sys,2);
            if norm(freqresp(sys,actFreq),2)>expNorm2
                expNorm2=actNorm2;
                expFreq=actFreq;
            end
            if (actFreq==0 &&expFreq==0)
                expNorm2=actNorm2;
            end
            expNorm = [expNorm1,expNorm2,expFreq,expNorm3];
            verification(testCase, actNorm, expNorm);
        end
        function testMIMObench(testCase)
            load('cdplayer.mat');
            sysSparse=sss(A,B,C);
            sys=ss(full(A),full(B),full(C),zeros(2,2));
            actNorm1=norm(sysSparse);
            [actNorm2,actFreq]=norm(sysSparse,inf);
            sysSparse.D=zeros(size(sysSparse.D));
            actNorm3=norm(sysSparse,2);
            actNorm = [actNorm1,actNorm2,actFreq,actNorm3];
            expNorm1=norm(sys);
            [expNorm2,expFreq]=norm(sys,inf);
            sys.D=zeros(size(sys.D));
            expNorm3=norm(sys,2);
            if norm(freqresp(sys,actFreq),2)>expNorm2
                expNorm2=actNorm2;
                expFreq=actFreq;
            end
            if (actFreq==0 &&expFreq==0)
                expNorm2=actNorm2;
            end
            expNorm = [expNorm1,expNorm2,expFreq,expNorm3];
            verification(testCase, actNorm, expNorm);
        end
        function testMIMOrandom(testCase)
            n=35;
            nInputs=7;
            nOutputs=5;
            sys=rss(n);
            sys=ss(sys.A,rand(n,nInputs),rand(nOutputs,n),rand(nOutputs,nInputs));
            sysSparse=sss(sys);
            actNorm1=norm(sysSparse);
            [actNorm2,actFreq]=norm(sysSparse,inf);
            sysSparse.D=zeros(size(sysSparse.D));
            actNorm3=norm(sysSparse,2);
            actNorm = [actNorm1,actNorm2,actFreq,actNorm3];
            expNorm1=norm(sys);
            [expNorm2,expFreq]=norm(sys,inf);
            sys.D=zeros(size(sys.D));
            expNorm3=norm(sys,2);
            if norm(freqresp(sys,actFreq),2)>expNorm2
                expNorm2=actNorm2;
                expFreq=actFreq;
            end
            if (actFreq==0 &&expFreq==0)
                expNorm2=actNorm2;
            end
            expNorm = [expNorm1,expNorm2,expFreq,expNorm3];
            verification(testCase, actNorm, expNorm);
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution(1:4),  expSolution(1:4),'RelTol',1e-3,...
    'Difference between actual and expected exceeds relative tolerance');
end