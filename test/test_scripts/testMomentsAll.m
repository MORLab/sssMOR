classdef testMomentsAll < matlab.unittest.TestCase
% testMoments - testing of moments.m using all benchmarks in local folder
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
% Last Change:  31 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------ 
    properties
         sysCell;
         Path;
    end
 
    methods(TestMethodSetup)
        function getBenchmarks(testCase)
            % change path
            testCase.Path = pwd;
            
            %insert path of local benchmark folder
            pathBenchmarks='M:\Dateien\Hiwi\testBenchmarks';
            cd(pathBenchmarks);

            % load files
            files = dir('*.mat'); 
            testCase.sysCell=cell(length(files));

            for i=1:length(files)
                eval(['load ' files(i).name]);
                if exist('E','var')
                    testCase.sysCell{i}=sss(A,B,C,0,E);
                    clear A B C E;
                elseif exist('A','var')
                    testCase.sysCell{i}=sss(A,B,C);
                    clear A B C;
                elseif exist('M','var')
                    E=blkdiag(speye(size(M)), M);
                    A=[zeros(size(M)),speye(size(M)); -K, -D];
                    B=[zeros(size(M,1),1); B];
                    C=[C, zeros(1,size(M,1))];     
                    testCase.sysCell{i}=sss(A,B,C,0,E);
                    clear A B C E;
                end
            end

            % change path back
            cd(testCase.Path);
        end
    end
    
    methods(Test)
        function testMoments1(testCase) 
            for i=1:length(testCase.sysCell)  
                if testCase.sysCell{i}.isDescriptor==0
                    sys=testCase.sysCell{i};
                    A=testCase.sysCell{i}.A;
                    B=testCase.sysCell{i}.B;
                    C=testCase.sysCell{i}.C;
                    s0=5;

                    m = moments(sys, s0, 4);
                    actSolution={m};

                    expm=zeros(1,4);
                    expm(1)=C*(A-s0*eye(size(A)))^(-1)*B;
                    expm(2)=C*(A-s0*eye(size(A)))^(-2)*B;
                    expm(3)=C*(A-s0*eye(size(A)))^(-3)*B;
                    expm(4)=C*(A-s0*eye(size(A)))^(-4)*B;

                    expSolution={expm};
                    verification(testCase, actSolution, expSolution, m);
                end
            end
        end
        
        function testMoments2(testCase) 
            %s0: imag
            for i=1:length(testCase.sysCell)  
                if testCase.sysCell{i}.isDescriptor==0
                    sys=testCase.sysCell{i};  
                    A=testCase.sysCell{i}.A;
                    B=testCase.sysCell{i}.B;
                    C=testCase.sysCell{i}.C;
                    s0=5+7i;

                    m = moments(sys, s0, 4);
                    actSolution={m};

                    expm=zeros(1,4);
                    expm(1)=C*(A-s0*eye(size(A)))^(-1)*B;
                    expm(2)=C*(A-s0*eye(size(A)))^(-2)*B;
                    expm(3)=C*(A-s0*eye(size(A)))^(-3)*B;
                    expm(4)=C*(A-s0*eye(size(A)))^(-4)*B;

                    expSolution={expm};
                    verification(testCase, actSolution, expSolution, m);
                end
            end
        end
        
        function testMoments3(testCase) 
            %s0: Inf
            for i=1:length(testCase.sysCell)  
                if testCase.sysCell{i}.isDescriptor==0
                    sys=testCase.sysCell{i};    
                    A=testCase.sysCell{i}.A;
                    B=testCase.sysCell{i}.B;
                    C=testCase.sysCell{i}.C;
                    s0=Inf;

                    m = moments(sys, s0, 4);
                    actSolution={m};

                    %mi=c*A^i*b
                    expm=zeros(1,2);
                    expm(1)=0; %m(1)=sys.D?
                    expm(2)=C*B;
                    expm(3)=C*A*B;
                    expm(4)=C*A^2*B;

                    expSolution={expm};
                    verification(testCase, actSolution, expSolution, m);
                end
            end
        end
        
         function testMoments4(testCase) 
            %with E-matrix, s0 real
            for i=1:length(testCase.sysCell)  
                if testCase.sysCell{i}.isDescriptor==1
                    sys=testCase.sysCell{i};    
                    A=testCase.sysCell{i}.A;
                    B=testCase.sysCell{i}.B;
                    C=testCase.sysCell{i}.C;
                    E=testCase.sysCell{i}.E;
                    s0=50;

                    m = moments(sys, s0, 4);
                    actSolution={m};

                    expm=zeros(1,4);
                    expm(1)=C*(A-s0*E)^(-1)*B;
                    expm(2)=C*(A-s0*E)^(-1)*E*(A-s0*E)^(-1)*B;
                    expm(3)=C*(A-s0*E)^(-1)*E*(A-s0*E)^(-1)*E*(A-s0*E)^(-1)*B;
                    expm(4)=C*(A-s0*E)^(-1)*E*(A-s0*E)^(-1)*E*(A-s0*E)^(-1)*E*(A-s0*E)^(-1)*B;

                    expSolution={expm};
                    verification(testCase, actSolution, expSolution, m);
                end  
            end
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