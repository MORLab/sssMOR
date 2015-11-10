classdef testIssd < matlab.unittest.TestCase
    
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
        function testIssd1(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                
                size(sys.A)
                [actSolution]=issd(sys);
                if nnz(eig(full(sys.A+sys.A'),'vector')>=0)
                    expSolution=0;
                else
                    expSolution=1;
                end
                if sys.isDescriptor
                    if nnz(eig(full(sys.E))<=0)
                        expSolution=0;
                    end
                end

                verification (testCase, actSolution, expSolution);
                verifyInstanceOf(testCase, actSolution , 'double', 'Instances not matching');
                verifySize(testCase, actSolution, [1 1], 'Size not matching');
            end
        end
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.1,'AbsTol',0.000001, ...
               'Difference between actual and expected exceeds relative tolerance');
end