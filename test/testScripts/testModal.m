classdef testModal < sssTest
% testModal - testing of modalMor.m
%
% Description:
%   The function modalMor.m is tested (3 tests) on:
%    + comparing the eigenvalues of the reduced system to the solution of
%      modreal (only 'SM' possible).
%    + test systems: diagonal, building, LF10 (with E-matrix)
%
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> morlab@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
%               Lisa Jeschek
% Last Change:  07 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------
 
    methods(Test)
        function testModal1(testCase) 
            %Diagonal matrix
            A=zeros(11);
            A=diag(-10:-10:-100);
            A(11,11)=-30;
            B=(1:11)';
            C=B';

            Opts.type='SM';
            [sysr] = modalMor(sss(A,B,C,0), 6, Opts);
            actSolution={full(sort(eig(sysr)))};
            
            [expsysr,~]=modreal(ss(full(A),full(B),full(C),0),6);
            expSolution={full(sort(eig(expsysr)))};
                     
            verification(testCase, actSolution, expSolution, sysr);
        end
        function testModal2(testCase) 
            %without E-matrix
            load('building.mat');

            Opts.type='SM';
            [sysr] = modalMor(sss(A,B,C,0), 6, Opts);
            actEig=sort(eig(sysr));
            actSolution={full(real(actEig)), full(abs(imag(actEig)))};
            
            [expsysr,~]=modreal(ss(full(A),full(B),full(C),0),6);
            expEig=sort(eig(expsysr));
            expSolution={full(real(expEig)), full(abs(imag(expEig)))};
                 
            verification(testCase, actSolution, expSolution, sysr);
        end
        function descriptorSymmetric(testCase) 
            %with E-matrix, symmetric
            sys = sss('rail_1357');
            
            q = 10;
            Opts.type='SM';
            [sysr] = modalMor(sys, q, Opts);
            actSolution=full(sort(eig(sysr)));            
            expSolution=full(sort(eigs(sys,q,Opts.type)));
                 
            verification(testCase, actSolution, expSolution, sysr);
        end
        function descriptor(testCase) 
            %with E-matrix, not symmetric
            warning('off'), sys = loadSss('LF10'); warning('on');
            
            q = 10;
            Opts.type = 'LM';
            [sysr] = modalMor(sys, q, Opts);
            actSolution=full(cplxpair(eig(sysr)));            
            expSolution=full(cplxpair(eigs(sys,q,Opts.type)));
                 
            verification(testCase, actSolution, expSolution, sysr);
        end
        function loopThroughOptions(testCase) 
            %run modal mor for all loaded models with all options to check
            %that all combinations are valid
            
            % Define possible opts combinations
            AllOpts.type    = {'SM','LM',1};
            AllOpts.orth    = {true, false,'qr'};
            AllOpts.real    = {true, false,'real'};
            AllOpts.tol     = {1e-6};
            AllOpts.dominance = {0,'analyze','2q','3q'};
            AllOpts.lse     = {'sparse','full'};
            AllOpts.subspaceW = {'eigs','1by1'};
            
            [AllOpts,nCases] = generateAllOpts(AllOpts);
            
            h = waitbar(0,'modalMor: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOpts{kOpts};
                for i=1:length(testCase.sysCell)
                    %  test system
                    sys     = testCase.sysCell{i};
                    if strcmp(Opts.dominance,'2q') && strcmp(testCase.sysCell{i}.Name,'iss')
                        % Option dominance = '2q' produces a numerical error
                        % for the benchmark "iss", because only 2 of the 
                        % requested 8 eigenvalues converge.
                        
                        warning('off')
                        try
                            sysr    = modalMor(sys,4, Opts);
                        end
                        warning('on')
                    else
                        sysr    = modalMor(sys,4, Opts);
                    end
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)                    
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution, sysr)
       verifyEqual(testCase, actSolution, expSolution,'RelTol',0.2,'AbsTol',0.00000001,...
            'Difference between actual and expected exceeds relative tolerance');
       verifyEqual(testCase, rank(full(sysr.A)), length(sysr.B),...
            'Rank(Ar) is not full');
       verifyEqual(testCase, rank(full(sysr.E)), length(sysr.B),...
            'Rank(Er) is not full');
       verifyEqual(testCase, nnz(isinf(sysr.A)), 0, ...
            'Ar contains Inf');
       verifyEqual(testCase, nnz(isinf(sysr.E)), 0, ...
            'Er contains Inf');
       verifyEqual(testCase, nnz(isnan(sysr.A)), 0, ...
            'Ar contains Nan');
       verifyEqual(testCase, nnz(isnan(sysr.E)), 0, ...
            'Er contains Nan');
       verifyLessThan(testCase, max(real(eig(sysr))), 0, ...
            'Ar is not stable');
end