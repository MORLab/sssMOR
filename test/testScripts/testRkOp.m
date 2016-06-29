classdef testRkOp < sssTest
    
    methods (Test)  
        function testRkOp1(testCase)
            Opts.rk='output';
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae && sys.isSiso
                    [h,t]=impulse(ss(sys));
                    sOpt1 = rkOp(h,t);
                    sOpt2 = rkOp(sys);
                    rkOp(ss(sys),10,Opts);
                    verification (testCase, sOpt1, sOpt2);
                end
            end
        end
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 2e-2,'AbsTol',1e-6, ...
               'Difference between actual and expected exceeds relative tolerance');
end