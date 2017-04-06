classdef testRkIcop < sssTest
    
    methods (Test)  
        function testRkIcop1(testCase)
            warning('off','sssMOR:rkIcop:sysrUnstable')
            warning('off','sssMor:rkIcop:NotConverged')
            Opts.rk='output';
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    sOpt1 = rkOp(sys);
                    success=1;
                    [~,~,~,sOpt2] = rkIcop(sys,10,30);
                    verification (testCase, sOpt1, sOpt2);
                    if success==1
                       rkIcop(ssRed(sys.A,sys.B,sys.C,sys.D,sys.E),0.5,12,Opts);
                    end
                end
            end
            warning('on','sssMOR:rkIcop:sysrUnstable')
            warning('on','sssMor:rkIcop:NotConverged')
        end
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 1e-2,'AbsTol',1e-6, ...
               'Difference between actual and expected exceeds relative tolerance');
end