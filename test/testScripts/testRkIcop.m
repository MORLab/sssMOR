classdef testRkIcop < sssTest
    
    methods (Test)  
        % FOLLOWING TEST WAS REMOVED AS RKICOP IS NOT GUARANTEED TO
        % CONVERGE TO THE OPTIMAL POINT
%         function compareRkOp2RkIcop(testCase)
%             warning('off','sssMOR:rkIcop:sysrUnstable')
%             warning('off','sssMor:rkIcop:NotConverged')
%             Opts.rk='output';
%             for i=1:length(testCase.sysCell)
%                 sys=testCase.sysCell{i};
%                 if ~sys.isDae
%                     % increase accuracy
%                     sOpt1       = rkOp(sys);
%                     
%                     Opts.tol    = 1e-5;
%                     Opts.maxiter= 1e2; 
%                     s0          = sOpt1; %put initial guess at optimum
%                     q           = 30;
%                     q           = ceil(q/(sys.p*sys.m));
%                     
%                     [~,~,~,sOpt2] = rkIcop(sys,s0,q,Opts);
%                     verification (testCase, sOpt2, sOpt1);
%                     
%                     rkIcop(ssRed(sys),s0,q,Opts);
%                 end
%             end
%             warning('on','sssMOR:rkIcop:sysrUnstable')
%             warning('on','sssMor:rkIcop:NotConverged')
%         end

         function runIcopWithoutVerification(testCase)
            warning('off','sssMOR:rkIcop:sysrUnstable')
            warning('off','sssMor:rkIcop:NotConverged')
            
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    % increase accuracy
                    sOpt1       = rkOp(sys);
                    
                    s0          = sOpt1; %put initial guess at optimum
                    q           = 30; q     = ceil(q/(sys.p*sys.m));
                    
                    [~,~,~,sOpt2] = rkIcop(sys,s0,q);                    
                    rkIcop(ssRed(sys),s0,q);
                end
            end
            warning('on','sssMOR:rkIcop:sysrUnstable')
            warning('on','sssMor:rkIcop:NotConverged')
        end
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 1e-1, ...
               'Difference between actual and expected exceeds relative tolerance');
end