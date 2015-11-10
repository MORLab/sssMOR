classdef testImpulse < matlab.unittest.TestCase
    
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
        function testImpulse1(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                t=1:0.1:5;
                
                [actH]=impulse(sys,t); 
                [expH]=impulse(ss(sys),t);
                         
                
                actSolution={actH};
                expSolution={expH};
                
                verification (testCase, actSolution, expSolution);
                verifyInstanceOf(testCase, actT , 'double', 'Instances not matching');
                verifyInstanceOf(testCase, actH , 'double', 'Instances not matching');
                verifySize(testCase, actH, size(expH), 'Size not matching');
                verifySize(testCase, actT, size(t), 'Size not matching');
            end
        end
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.1,'AbsTol',0.000001, ...
               'Difference between actual and expected exceeds relative tolerance');
end