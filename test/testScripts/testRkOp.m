classdef testRkOp < sssTest
    
    methods (Test)  
        function testRkOp1(testCase)
            Opts.rk='output';
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    [g,t]=impulse(ss(sys),0:0.01:decayTime(sys));
                    sOpt1 = rkOp(g,t);
                    sOpt2 = rkOp(sys);
                    rkOp(ssRed('userDefined',[],sys.A,sys.B,sys.C, ...
                         sys.D,sys.E),10,Opts);
                    verification (testCase, sOpt1, sOpt2);
                end
            end
        end
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 85e-2,'AbsTol',1e-6, ...
               'Difference between actual and expected exceeds relative tolerance');
end