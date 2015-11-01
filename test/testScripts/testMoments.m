classdef testMoments < matlab.unittest.TestCase
% testMoments - testing of moments.m
%
% Description:
%   The function moments.m is tested (4 tests) on:
%    + moments: mi=c*(A-s0*I)^-(i+1)*b i=0,1,...
%    + markov: mi=c*A^i*b, i=0,1,..
%    + test systems: build, beam, random, LF10 (with E-matrix)
%    + Neither Inf nor NaN in m
%    + s0: real, imag, Inf
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
% Last Change:  28 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------ 
    properties
    end
 
    
    methods(Test)
        function testMoments1(testCase) 
            %s0: real
            load('build.mat');      
            s0=5;
            
            m = moments(sss(A,B,C), s0, 4);
            actSolution={m};
            
            expm=zeros(1,1,4);
            expm(1)=C*(A-s0*eye(size(A)))^(-1)*B;
            expm(2)=C*(A-s0*eye(size(A)))^(-2)*B;
            expm(3)=C*(A-s0*eye(size(A)))^(-3)*B;
            expm(4)=C*(A-s0*eye(size(A)))^(-4)*B;
            
            expSolution={expm};
            verification(testCase, actSolution, expSolution, m);
        end
        
        function testMoments2(testCase) 
            %s0: imag
            load('beam.mat');      
            s0=5+7i;
            
            m = moments(sss(A,B,C), s0, 4);
            actSolution={m};
            
            expm=zeros(1,1,4);
            expm(1)=C*(A-s0*eye(size(A)))^(-1)*B;
            expm(2)=C*(A-s0*eye(size(A)))^(-2)*B;
            expm(3)=C*(A-s0*eye(size(A)))^(-3)*B;
            expm(4)=C*(A-s0*eye(size(A)))^(-4)*B;
            
            expSolution={expm};
            verification(testCase, actSolution, expSolution, m);
        end
        
        function testMoments3(testCase) 
            %s0: Inf
            load('random.mat');      
            s0=Inf;
            
            m = moments(sss(A,B,C,0), s0, 4);
            actSolution={m};
            
            %mi=c*A^i*b
            expm=zeros(1,1,4);
            expm(1)=0; %m(1)=sys.D?
            expm(2)=C*B;
            expm(3)=C*A*B;
            expm(4)=C*A^2*B;

            expSolution={expm};
            verification(testCase, actSolution, expSolution, m);
        end
        
        function testMoments4(testCase) 
            %with E-matrix, s0 real
            load('LF10.mat');
            E=blkdiag(speye(size(M)), M);
            A=[zeros(size(M)),speye(size(M)); -K, -D];
            B=[zeros(size(M,1),1); B];
            C=[C, zeros(1,size(M,1))];     
            s0=50;
            
            m = moments(sss(A,B,C,0,E), s0, 4);
            actSolution={m};
            
            expm=zeros(1,1,4);
            expm(1)=C*(A-s0*E)^(-1)*B;
            expm(2)=C*(A-s0*E)^(-1)*E*(A-s0*E)^(-1)*B;
            expm(3)=C*(A-s0*E)^(-1)*E*(A-s0*E)^(-1)*E*(A-s0*E)^(-1)*B;
            expm(4)=C*(A-s0*E)^(-1)*E*(A-s0*E)^(-1)*E*(A-s0*E)^(-1)*E*(A-s0*E)^(-1)*B;
   
            expSolution={expm};
            verification(testCase, actSolution, expSolution, m);
        end
    end
end

function [] = verification(testCase, actSolution, expSolution, m)
       verifyEqual(testCase, actSolution, expSolution,'RelTol',0.1,...
            'Difference between actual and expected exceeds relative tolerance');
       verifyEqual(testCase, nnz(isinf(m)), 0, ...
            'm contains Inf');
       verifyEqual(testCase, nnz(isnan(m)), 0, ...
            'm contains Nan');
end