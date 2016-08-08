classdef testIsrk < sssTest
    methods(Test)
        function testIsrk1(testCase) 
            Opts=struct('tol',0.05,'maxiter',100,'stopCrit','combAny');
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae && sys.isSiso
                    
                    s0=[0,0,100,1+5i,1-5i,14-0.2i,14+0.2i, Inf, Inf];
                    [sysr_old, V, W, s0_old, Rt, B_, Ssylv, Rsylv, kIter, s0_Traj, RtTraj, flag] = isrk(sys, s0, Opts);
                    
                    if flag==0
                        expV=arnoldi(sys.E,sys.A,sys.B,s0_old);
                        Q = lyap(sys.A',sys.C'*sys.C,[],sys.E');
                        expW=Q*sys.E*expV;
                        
                        
                        sysr=sss(expW'*sys.A*expV, expW'*sys.B, sys.C*expV, sys.D, expW'*sys.E*expV);
                        s0=-eig(sysr)';
                        
                        stopCrit = norm((s0-s0_old)./s0, 1)/sysr.n;
                        stopCrit = [stopCrit, inf]; 
                        if all(real(eig(sysr))<0) && all(real(eig(sysr_old))<0)
                                stopCrit(2) = norm(sysr-sysr_old)/norm(sysr);
                        end
                        stopCrit=min(stopCrit);

                        
                        actM = moments(sysr_old,s0_old, 1); expM = moments(sys,s0_old, 1);
                        expSolution={actM};
                        actSolution={expM};

                        verification(testCase, actSolution, expSolution, sysr_old, s0_old, s0_Traj);
                        verifyEqual(testCase, isstable(sysr), 1, 'Sysr is not stable');
                        verifyLessThan(testCase, stopCrit, Opts.tol);
                    
                    else
                        sys.Name
                        warning(['Isrk did not converge within ', num2str(Opts.maxiter),' steps.']);    
                    end
                end
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution, sysr, s0, s0_traj)
       verifyEqual(testCase, actSolution, expSolution,'RelTol',0.3,'AbsTol',0.00000001,...
            'Difference between actual and expected exceeds relative tolerance');
       verifyLessThanOrEqual(testCase, max(imag(sysr.A)), 0, ...
            'Ar is not purely real'); 
       verifyLessThanOrEqual(testCase, max(imag(sysr.E)), 0, ...
            'Er is not purely real'); 
       verifyEqual(testCase, rank(full(sysr.A)), length(s0),...
            'Rank(Ar) is not full');
       verifyEqual(testCase, rank(full(sysr.E)), length(s0),...
            'Rank(Er) is not full');
       verifyEqual(testCase, nnz(isinf(sysr.A)), 0, ...
            'Ar contains Inf');
       verifyEqual(testCase, nnz(isinf(sysr.E)), 0, ...
            'Er contains Inf');
       verifyEqual(testCase, nnz(isnan(sysr.A)), 0, ...
            'Ar contains Nan');
       verifyEqual(testCase, nnz(isnan(sysr.E)), 0, ...
            'Er contains Nan');
       verifyLessThanOrEqual(testCase, size(s0_traj,1),301, ...
            'k > kmax');
end