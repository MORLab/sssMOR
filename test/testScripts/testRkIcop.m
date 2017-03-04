classdef testRkIcop < sssTest
    
    methods (Test)  
        function testRkIcop1(testCase)
            Opts.rk='output';
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    sOpt1 = rkOp(sys);
                    success=1;
                    try
                        [~,~,~,sOpt2] = rkIcop(sys,10,30);
                    catch err
                        if strcmp(err.identifier,'sssMor:rkIcopNotConverged')
                            success=0;
                        else
                            error(err.message);
                        end
                    end
                    if success==1
                        try
                            rkIcop(ssRed(sys.A,sys.B,sys.C,sys.D,sys.E),0.5,12,Opts);
                        catch err
                            if strcmp(err.identifier,'sssMor:rkIcopNotConverged')
                                success=0;
                            else
                                error(err.message);
                            end
                        end
                        if success==1
                            verification (testCase, sOpt1, sOpt2);
                        end
                    end
                end
            end
        end
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 1e-2,'AbsTol',1e-6, ...
               'Difference between actual and expected exceeds relative tolerance');
end