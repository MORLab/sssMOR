classdef testConnect < matlab.unittest.TestCase
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
    methods(Test)  
        function testConnect1(testCase)
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
            if sys1.m>5 || sys2.m>5 || sys1.p>5 || sys2.p>5
                error('test is only for systems with sys.p<6 and sys.m<6');
            end
            inputVector1={'u11','u12','u13','u14','u15'};
            outputVector1={'y11','y12','y13','y14','y15'};
            inputVector2={'u21','u22','u23','u24','u25'};
            outputVector2={'y21','y22','y23','y24','y25'};
            sys1.InputName=inputVector1(1:sys1.m)';
            sys1.OutputName=outputVector1(1:sys1.p)';
            sys2.InputName=inputVector2(1:sys2.m)';
            sys2.OutputName=outputVector2(1:sys2.p)';

            [actSys]=connect(sys1,sys2,[inputVector1(1:sys1.m),inputVector2(1:sys2.m)],[outputVector1(1:sys1.p),outputVector2(1:sys2.p)]);
            [expSys]=connect(ss(sys1),ss(sys2),[inputVector1(1:sys1.m),inputVector2(1:sys2.m)],[outputVector1(1:sys1.p),outputVector2(1:sys2.p)]);

            actSolution={full(actSys.E\actSys.A), full(actSys.B), full(actSys.C), full(actSys.D)};
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

