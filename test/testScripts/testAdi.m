classdef testAdi < matlab.unittest.TestCase
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
% Last Change:  05 Sep 2015
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
        function testAdi1(testCase) %tbr
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae && sys.n>100
                    q=50;
                    opts.adi='adi';
                    tbr(sys,q,opts);
                    actLambda=(svd(full(sys.ObsGramChol'*sys.ConGramChol)));
                    
                    actSolution={actLambda(1:5)/max(actLambda)};
                    
                    opts.adi=0;
                    sys.ConGramChol=[];
                    sys.ObsGramChol=[];
                    tbr(sys,q,opts);
                    expLambda=(svd(full(sys.ConGramChol*sys.ObsGramChol')));
                    expSolution={expLambda(1:5)/max(expLambda)};

                    verifyEqual(testCase, actSolution, expSolution,'RelTol',0.3,...
                        'Difference between actual and expected exceeds relative tolerance');
                end
            end
        end
%         function testAdi2(testCase) %norm
%             for i=1:length(testCase.sysCell)
%                 sys=testCase.sysCell{i};
%                 if ~sys.isDae && sys.n>100
%                     q=50;
%                     opts.adi='adi';
%                     norm(sys,q,opts);
%                     actLambda=(svd(full(sys.ObsGramChol'*sys.ConGramChol)));
%                     
%                     actSolution={actLambda(1:5)/max(actLambda)};
%                     
%                     opts.adi=0;
%                     sys.ConGramChol=[];
%                     sys.ObsGramChol=[];
%                     norm(sys,q,opts);
%                     expLambda=(svd(full(sys.ConGramChol*sys.ObsGramChol')));
%                     expSolution={expLambda(1:5)/max(expLambda)};
% 
%                     verifyEqual(testCase, actSolution, expSolution,'RelTol',0.3,...
%                         'Difference between actual and expected exceeds relative tolerance');
%                 end
%             end
%         end
    end
end