classdef testIsstable < matlab.unittest.TestCase
    
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
    
    %Test functions
    methods (Test)  
        function testIsstable1(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                
                [actB] = isstable(sys);
                [expB]=isstable(ss(sys));
                
                actSolution={logical(actB)};
                expSolution={expB};
                
                verification (testCase, actSolution, expSolution);
            end
        end
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.1,'AbsTol',0.000001, ...
               'Difference between actual and expected exceeds relative tolerance');
end