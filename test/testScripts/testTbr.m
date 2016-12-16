classdef testTbr < sssTest
% testtbr - testing of tbr.m
%
% Description:
%   The function tbr.m is tested (5 tests) on:
%    + q: 5, 10, 15, 20, 25
%    + test systems: building, beam, fom, random, LF10 (with E-matrix).
%    + comparing sysr with directly calculated solution 
%    + W'*V identity matrix
%    + Ar, Er purely real
%    + rank(Ar), rank(Er) full
%    + Neither Inf nor NaN in Ar, Er
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
%               Lisa Jeschek
% Last Change:  11 Feb 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------ 
    
    methods(Test)
        function testTbr1(testCase) %q=5
            sys = loadSss('building.mat');      
            q=5;
            
            [sysr, V, W, hsv, S, R] = tbr(sys,q);            
            verification(testCase,sys, sysr, V, W, hsv, S, R);
        end
        
        function testBalancedRealizationExplicit(testCase)
            sys = loadSss('building.mat');      
            
            [sysr, V, W, hsv, S, R] = tbr(sys,sys.n);        
            verification(testCase,sys, sysr, V, W, hsv, S, R);
        end
        function testBalancedRealizationImplicit(testCase)
            sys = loadSss('building.mat');    
            sys = rk(sys,zeros(1,10),zeros(1,10)); %create model with E~=I
            
            [sysr, V, W, hsv, S, R] = tbr(sys,sys.n);        
            verification(testCase,sys, sysr, V, W, hsv, S, R);
        end
        
        function testTbr2(testCase) %q=15
            sys = loadSss('beam.mat');
            q=15;
            
            [sysr, V, W, hsv, S, R] = tbr(sys,q);
            
            verification(testCase,sys, sysr, V, W, hsv, S, R);

        end        
        function testTbr3(testCase) %q=25
            sys = loadSss('fom.mat');
            q=25;
            Opts.type='tbr';
            
            [sysr, V, W, hsv, S, R] = tbr(sys,q, Opts);
            verification(testCase,sys, sysr, V, W, hsv, S, R);
        end
        function testTbr5(testCase) %q=10 (with E-matrix)
            load('LF10.mat');
            q=10;
            
            E=blkdiag(speye(size(M)), M);
            A=[zeros(size(M)),speye(size(M)); -K, -D];
            B=[zeros(size(M,1),1); B];
            C=[C, zeros(1,size(M,1))];
            
            sys = sss(A,B,C,[],E);
            
            [sysr, V, W, hsv, S, R] = tbr(sys,q);
            verification(testCase,sys, sysr, V, W, hsv, S, R);
        end
        function testTbr6(testCase) %adi
            sys = loadSss('rail_1357');
            q=50;
            opts.type='adi';
            [sysr,V,W,hsv,S,R]=tbr(sys,q,opts);
            verification(testCase,sys, sysr, V, W, hsv, S, R);
        end
        function testMatchDcGain(testCase)
            warning('on','tbr:rcond');
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    Opts.type='matchDcGain';
                    sysr=tbr(sys,10,Opts);
                    w = warning('query','last');
                    if isempty(w) || ~strcmp(w.identifier,'tbr:rcond') %A22 not close to singular
                        actSolution= freqresp(sysr,0);
                        expSolution= freqresp(sys,0);
                        verifyEqual(testCase, actSolution , expSolution, 'AbsTol', 1e-3, 'RelTol', 1e-3);
                    end
                end
            end
        end
    end
end

function [] = verification(testCase, sys, sysr,V,W,hsvs,S,R)
        X = S*S'; Y = R*R';  
        tol = 1e-3;
       
       % Solution of Lyapunov equations
       verifyLessThan(testCase, ...
           norm(sys.A*(X)*sys.E' + sys.E*(X)*sys.A' + sys.B*sys.B'),...
           tol,'C-Lyapunov Equation not satisfied');           
       verifyLessThan(testCase, ...
            norm(sys.A'*(Y)*sys.E + sys.E'*(Y)*sys.A + sys.C'*sys.C),...
            tol,'O-Lyapunov Equation not satisfied');  
        
        % Gramians
        if sys.n < 500
           P = gram(ss(sys),'c'); Q = gram(ss(sys),'o');
           verifyLessThan(testCase, norm(X-P), tol,'Controllability Gramian')
           verifyLessThan(testCase, norm(sys.E'*Y*sys.E-Q), tol,'Observability Gramian')
        end
       
       % Balanced ROM
       Pr = gram(sysr,'c'); Qr = gram(sysr,'o');
       verifyLessThan(testCase, norm(diag(diag(Pr))-Pr),tol,'ROM not balanced');
       verifyLessThan(testCase, norm(diag(diag(Qr))-Qr),tol,'ROM not balanced');
       
       % ROM HVS
       hsvr = sqrt(eig(Pr*Qr));
       verifyLessThan(testCase, sort(hsvs(1:sysr.n))-sort(hsvr), tol, 'HSV in ROM don''t match FOM')
       
       % Stable ROM
       verifyTrue(testCase,isstable(sysr),'ROM unstable');


       % Biorthogonal W,V wrt E
       if sysr.n<sys.n %approximation
        verifyLessThan(testCase, norm((W'*sys.E*V)-eye(sysr.n)),tol,...
               'W''*E*V not identity matrix');
       else %balanced realization
           verifyLessThan(testCase, norm((W'*V)-eye(sysr.n)),tol,...
               'W''*V not identity matrix');
       end
           
end