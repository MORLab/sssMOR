classdef testConnectSss < matlab.unittest.TestCase
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
    methods(Test)  
        function testConnectSss1(testCase)
            if length(testCase.sysCell)>1
                for i=1:2:length(testCase.sysCell)-1
                    sys1=testCase.sysCell{i};
                    sys2=testCase.sysCell{i+1};
                    sysConnect(testCase,sys1,sys2)
                end
            else 
                sys1=loadSss('build.mat');
                sys2=loadSss('random.mat');
                sysConnect(testCase,sys1,sys2);
            end
        end
    end
end

function []=sysConnect(testCase,sys1,sys2)
            sys=append(sys1,sys2);
            K=rand(sys.p, sys.m);

            [actSys]=connectSss(sys,K);
            [expSys]=feedback(ss(sys),-K);

            actSolution={full(actSys.A), full(actSys.B), full(actSys.C), full(actSys.D)};
            expSolution={expSys.A, expSys.B, expSys.C, expSys.D};

            verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.1,'AbsTol',0.000001, ...
               'Difference between actual and expected exceeds relative tolerance');
            verifyInstanceOf(testCase, full(actSys.A) , 'double', 'Instances not matching');
            verifyInstanceOf(testCase, full(actSys.B) , 'double', 'Instances not matching');
            verifyInstanceOf(testCase, full(actSys.C) , 'double', 'Instances not matching');
            verifyInstanceOf(testCase, full(actSys.D) , 'double', 'Instances not matching');
            verifyInstanceOf(testCase, full(actSys.E) , 'double', 'Instances not matching');
            verifySize(testCase, actSys.A, size(expSys.A), 'Size not matching');
            verifySize(testCase, actSys.B, size(expSys.B), 'Size not matching');
            verifySize(testCase, actSys.C, size(expSys.C), 'Size not matching');
            verifySize(testCase, actSys.D, size(expSys.D), 'Size not matching');
            verifySize(testCase, actSys.E, size(expSys.A), 'Size not matching');
end

