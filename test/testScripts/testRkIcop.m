classdef testRkIcop < sssTest
    
    methods (Test)  
        function testRkIcop1(testCase)
            Opts.rk='output';
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    sOpt1 = rkOp(sys);
                    [~,~,~,sOpt2] = rkIcop(sys,10,30);
                    rkIcop(ss(sys),0.5,12,Opts);
                    verification (testCase, sOpt1, sOpt2);
                end
            end
        end
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 1e-2,'AbsTol',1e-6, ...
               'Difference between actual and expected exceeds relative tolerance');
end