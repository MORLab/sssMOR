classdef testDiag < matlab.unittest.TestCase
    
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
    
    %Test functions
    methods (Test)  
        function testDiag1(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                
                [actSys]=diag(sys);
                [expSys]=canon(ss(sys),'modal');
                
                actSolution={trace(full(actSys.A)),sort(real(eig(full(actSys.A)))), sort(imag(eig(full(actSys.A))))};
                expSolution={trace(expSys.A),sort(real(eig(expSys.A))), sort(imag(eig(expSys.A)))};
                
                verification (testCase, actSolution, expSolution);
                verifyInstanceOf(testCase, full(actSys.A) , 'double', 'Instances not matching');
                verifyInstanceOf(testCase, full(actSys.B) , 'double', 'Instances not matching');
                verifyInstanceOf(testCase, full(actSys.C) , 'double', 'Instances not matching');
                verifyInstanceOf(testCase, full(actSys.E) , 'double', 'Instances not matching');
                verifySize(testCase, actSys.A, size(sys.A), 'Size not matching');
                verifySize(testCase, actSys.B, size(sys.B), 'Size not matching');
                verifySize(testCase, actSys.C, size(sys.C), 'Size not matching');
                verifySize(testCase, actSys.E, size(sys.E), 'Size not matching');
                verifyLessThan(testCase, norm(full(actSys.E)-eye(size(actSys.E))),0.01, 'E not identity');
            end
        end
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.1,'AbsTol',0.000001, ...
               'Difference between actual and expected exceeds relative tolerance');
end