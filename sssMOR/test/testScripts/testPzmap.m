classdef testPzmap < matlab.unittest.TestCase
    % testPzmap - testing of pzmap.m
    %
    % Description:
    %   The function pzmap.m is tested (3 tests) on:
    %    + Norm of a SISO benchmark system.
    %    + Norm of a SISO random system.
    %    + Norm of a MISO random system.
    %    + Norm of a SIMO random system.
    %    + Norm of MIMO benchmark system.
    %    + Norm of a MIMO random system.
    %    + Verifies for every case the following inputs/outputs:
    %    pzmap(sys),[p,z]=pzmap(sys)
    %
    %
    % ------------------------------------------------------------------
    %   This file is part of sssMOR, a Sparse State Space, Model Order
    %   Reduction and System Analysis Toolbox developed at the Institute
    %   of Automatic Control, Technische Universitaet Muenchen.
    %   For updates and further information please visit www.rt.mw.tum.de
    %   For any suggestions, submission and/or bug reports, mail us at
    %                     -> sssMOR@rt.mw.tum.de <-
    % ------------------------------------------------------------------
    % Authors:      Alessandro Castagnotto, Maria Cruz Varona
    %               Jorge Luiz Moreira Silva
    % Last Change:  05 Nov 2015
    % Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
    % ------------------------------------------------------------------
    methods(Test)
        function testSISObench(testCase)
            load('build.mat');
            sysSparse=sss(A,B,C);
            sys=ss(A,B,C,zeros(1,1));
            [actP,actZ]=pzmap(sysSparse);
            actPZ={sort(actP),sort(actZ)};
            pzmap(sysSparse);
            [expP,expZ]=pzmap(sys);
            expPZ={sort(expP),sort(expZ)};
            verification(testCase, actPZ, expPZ);
            close all;
        end
        function testSISOrandom(testCase)
            sys=rss(35);
            sysSparse=sss(sys);
            [actP,actZ]=pzmap(sysSparse);
            actPZ={sort(actP),sort(actZ)};
            pzmap(sysSparse);
            [expP,expZ]=pzmap(sys);
            expPZ={sort(expP),sort(expZ)};
            verification(testCase, actPZ, expPZ);
            close all;
        end
        function testMISOrandom(testCase)
            n=35;
            nInputs=5;
            sys=rss(n);
            sys=ss(sys.A,rand(n,nInputs),sys.C,rand(1,nInputs));
            sysSparse=sss(sys);
            [actP,actZ]=pzmap(sysSparse);
            actPZ={sort(actP),sort(actZ)};
            pzmap(sysSparse);
            [expP,expZ]=pzmap(sys);
            expPZ={sort(expP),sort(expZ)};
            verification(testCase, actPZ, expPZ);
            close all;
        end
        function testSIMOrandom(testCase)
            n=35;
            nOutputs=5;
            sys=rss(n);
            sys=ss(sys.A,sys.B,rand(nOutputs,n),rand(nOutputs,1));
            sysSparse=sss(sys);
            [actP,actZ]=pzmap(sysSparse);
            actPZ={sort(actP),sort(actZ)};
            pzmap(sysSparse);
            [expP,expZ]=pzmap(sys);
            expPZ={sort(expP),sort(expZ)};
            verification(testCase, actPZ, expPZ);
            close all;
        end
        function testMIMObench(testCase)
            load('cdplayer.mat');
            sysSparse=sss(A,B,C);
            sys=ss(full(A),full(B),full(C),zeros(2,2));
            [actP,actZ]=pzmap(sysSparse);
            actPZ={sort(actP),sort(actZ)};
            pzmap(sysSparse);
            [expP,expZ]=pzmap(sys);
            expPZ={sort(expP),sort(expZ)};
            verification(testCase, actPZ, expPZ);
            close all;
        end
        function testMIMOrandom(testCase)
            n=35;
            nInputs=7;
            nOutputs=5;
            sys=rss(n);
            sys=ss(sys.A,rand(n,nInputs),rand(nOutputs,n),rand(nOutputs,nInputs));
            sysSparse=sss(sys);
            [actP,actZ]=pzmap(sysSparse);
            actPZ={sort(actP),sort(actZ)};
            pzmap(sysSparse);
            [expP,expZ]=pzmap(sys);
            expPZ={sort(expP),sort(expZ)};
            verification(testCase, actPZ, expPZ);
            close all;
        end

    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',0.1e-6,...
    'Difference between actual and expected exceeds relative tolerance');
end