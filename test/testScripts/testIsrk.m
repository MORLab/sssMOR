classdef testIsrk < sssTest
    methods(Test)
        function testIsrk1(testCase) 
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae && sys.isSiso
                    Opts=struct('epsilon',0.05,'maxiter',50,'stopCrit','s0','real',true);
                    s0=[0,0,50,50,100,1+5i,1-5i,14-0.2i,14+0.2i, Inf, Inf];
                    [sysr, V, W, s0, Rt, B_, Ssylv, Rsylv, kIter, s0_Traj, RtTraj, flag] = isrk(sys, s0, Opts);
                    
                    if flag==0
                        expV=arnoldi(sys.E,sys.A,sys.B,sort(s0));
                        A = sys.E\sys.A;
                        Q = lyap(A',sys.C'*sys.C);
                        expW=Q*expV/(expV'*Q*expV);
                        
                        actM = moments(sysr,s0, 1); expM = moments(sys,s0, 1);

                        expSolution={expW'*A*expV, expW'*sys.B, sys.C*expV, expV, expW, actM};
                        actSolution={full(sysr.E\sysr.A), full(sysr.E\sysr.B), full(sysr.C),V, W, expM};

                        verification(testCase, actSolution, expSolution, sysr, s0, s0_Traj);
                        verifyEqual(testCase, isstable(sysr), 1, 'Sysr is not stable');
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