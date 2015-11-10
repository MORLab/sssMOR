classdef testIsstable < matlab.unittest.TestCase
    
    properties
         sysCell;
    end
 
    methods(TestMethodSetup)
        function getBenchmarks(testCase)
            load('benchmarksSysCell.mat');
            testCase.sysCell=benchmarksSysCell;
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