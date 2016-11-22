classdef testCirka < sssTest
% testIrka - testing of cirka.m
%
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Last Change:  22 Nov 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------ 
 
    methods(Test)
        function generalFunctionality(testCase) 
            warning('off','sssMOR:irka:maxiter')
            n = 2;
            s0 = zeros(1,n);
            
            for i=1:length(testCase.sysCell)
                %  test system
                sys=testCase.sysCell{i};  
                sys= sys(1,1);
                
                [sysr, V, W, s0, kIrka, sysm, relH2err] = cirka(sys, s0);
                
                verifyClass(testCase,sysr,'ssRed')
                verifySize(testCase,V,[sys.n,n])
                verifySize(testCase,W,[sys.n,n])
%                 verifyEqual(testCase,sysr,projectiveMor(sys,V,W),'RelTol',1e-3)
                verifyTrue(testCase,isstable(sysm))
                verifyGreaterThanOrEqual(testCase,size(sysm.A,1),n);
                verifyEqual(testCase,norm(sysm-sysr)/norm(sysm),relH2err)           
            end
            warning('on','sssMOR:irka:maxiter')
        end
        function convergenceCriteria(testCase) 
            warning('off','sssMOR:irka:maxiter')
            n = 2;
            s0 = zeros(1,n);
            
            crit = {'s0','sysr','sysm'};
            for j = 1:length(crit)
                Opts.stopcrit = crit{j};
                Opts.tol      = 1e-5;   
                Opts.maxiter  = 10;
                for i=1:length(testCase.sysCell)
                    %  test system
                    sys=testCase.sysCell{i};  
                    sys= sys(1,1);

                    [sysr, ~, ~, s0opt, kIrka, ~, ~] = cirka(sys, s0,Opts);
                    if length(kIrka)< Opts.maxiter
                        if strcmp(crit{j},'s0')
                            verifyEqual(testCase,-s0opt',eig(sysr),'RelTol',Opts.tol)
                            verifyTrue(testCase,isH2opt(sys,sysr,s0opt,struct('tol',Opts.tol)))
                        end
                    end
                end
            end
            warning('on','sssMOR:irka:maxiter')
        end
    end
end