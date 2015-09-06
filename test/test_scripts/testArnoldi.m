classdef testArnoldi < matlab.unittest.TestCase
% testArnoldi - testing of arnoldi.m
%
% Description:
%   The function arnoldi.m is tested (5 tests) on:
%    + comparing the columns of V to a directly calculated solution 
%    + AV=A*V, EV=E*V
%    + V'*V orthonormal
%    + V purely real
%    + rank(V) full
%    + Neither Inf nor NaN in V
%    + s0: purely real, purely imaginary, zero, Inf (Markov-parameter)
%    + test systems: build, beam, fom, random, LF10 (with E-matrix)
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
%               Lisa Jeschek
% Last Change:  04 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------
    properties
    end
    
    methods(Test)        
        function testArnoldi1(testCase) 
            %s0=Inf, real, imag (without E-matrix)
            load('build.mat');

            [V] = arnoldi(speye(size(A)),A,B,[Inf, 50, 100, 200, 300, 1-1i, 1+1i]);
            actSolution={V};
            
            temp=full(real(B));
            expV1=temp/norm(temp);      
            
            temp=(A-50*eye(size(A)))\B;
            expV2=temp/norm(temp);
            
            temp=(A-100*eye(size(A)))\B;
            expV3=temp/norm(temp);
            
            temp=(A-200*eye(size(A)))\B;
            expV4=temp/norm(temp);
            
            temp=(A-300*eye(size(A)))\B;
            expV5=temp/norm(temp);
            
            temp=real((A-(1-1i)*eye(size(A)))\B);
            expV6=temp/norm(temp);
            
            temp=imag((A-(1-1i)*eye(size(A)))\B);
            expV7=temp/norm(temp);
     
            [expV,~]=qr([expV1, expV2, expV3, expV4, expV5, expV6, expV7]);
            expV=[-expV(:,1),expV(:,2), expV(:,3), -expV(:,4), -expV(:,5), -expV(:,6), expV(:,7)];
            expSolution={expV};

            verification(testCase, actSolution,expSolution,V)
        end
        
        function testArnoldi2(testCase) 
            %s0: only imag (multiple value)
            load('beam.mat');

            [V] = arnoldi(speye(size(A)),A,B,[1-1i, 1-1i, 1-1i, 1+1i, 1+1i, 1+1i]);
            actSolution={V};
            
            temp=real((A-(1-1i)*eye(size(A)))\B);
            expV1=temp/norm(temp);
 
            temp=real((A-(1-1i)*eye(size(A)))\expV1);
            expV2=temp/norm(temp);
            
            temp=real((A-(1-1i)*eye(size(A)))\expV2);
            expV3=temp/norm(temp);
            
            temp=imag((A-(1-1i)*eye(size(A)))\B);
            expV4=temp/norm(temp);
            
            temp=imag((A-(1-1i)*eye(size(A)))\expV1);
            expV5=temp/norm(temp);
            
            temp=imag((A-(1-1i)*eye(size(A)))\expV2);
            expV6=temp/norm(temp);
            
            [expV,~]=qr([expV1, expV2, expV3, expV4, expV5, expV6]);
            expV=[expV(:,1), -expV(:,2), expV(:,3), -expV(:,4),expV(:,5), -expV(:,6)];
            expSolution={expV};

            verification(testCase, actSolution,expSolution,V)
        end
        function testArnoldi3(testCase) 
            %s0: only real (multiple value)
            load('random.mat');

            [V] = arnoldi(speye(size(A)),A,B,[50,50,50,50]);
            actSolution={V};
     
            temp=(A-50*eye(size(A)))\B;
            expV1=temp/norm(temp);
            
            temp=(A-50*eye(size(A)))\expV1;
            expV2=temp/norm(temp);
            
            temp=(A-50*eye(size(A)))\expV2;
            expV3=temp/norm(temp);
            
            temp=(A-50*eye(size(A)))\expV3;
            expV4=temp/norm(temp);
            
            [expV,~]=qr([expV1, expV2, expV3, expV4]);
            expV=[-expV(:,1), -expV(:,2), expV(:,3), expV(:,4)];
            expSolution={expV};

            verification(testCase, actSolution,expSolution,V)
        end
        
        function testArnoldi4(testCase) 
            %s0: only Inf (multiple value)
            load('fom.mat');

            [V] = arnoldi(speye(size(A)),A,B,[Inf, Inf, Inf, Inf]);
            actSolution={V};
     
            temp=full(B);
            expV1=temp/norm(temp);
            
            temp=full(A*B);
            expV2=temp/norm(temp);

            temp=full(A^2*B);
            expV3=temp/norm(temp);
            
            temp=full(A^3*B);
            expV4=temp/norm(temp);
            
            [expV,~]=qr([expV1, expV2, expV3, expV4]);
            expV=[-expV(:,1), -expV(:,2), -expV(:,3), expV(:,4)];
            expSolution={expV};

            verification(testCase, actSolution,expSolution,V)
        end
        function testArnoldi5(testCase) 
            %s0=Inf, real, imag, zero with E-matrix
            load('LF10.mat');
            E=blkdiag(speye(size(M)), M);
            A=[zeros(size(M)),speye(size(M)); -K, -D];
            B=[zeros(size(M,1),1); B];

            [V] = arnoldi(E,A,B,[Inf, 0, 100, 4+13i, 4-13i], @(x,y) (x'*y));
            actSolution={V};
             
            temp=full(real(E\B));
            expV1=temp/norm(temp);      
            
            temp=full((A-0*E)\B);
            expV2=temp/norm(temp);
            
            temp=full((A-100*E)\B);
            expV3=temp/norm(temp);
            
            temp=full(real((A-(4+13i)*E)\B));
            expV4=temp/norm(temp);
            
            temp=full(imag((A-(4+13i)*E)\B));
            expV5=temp/norm(temp);
     
            [expV,~]=qr([expV1, expV2, expV3, expV4, expV5]);
            expV=[-expV(:,1), expV(:,2), expV(:,3), expV(:,4), expV(:,5)];
            expSolution={expV};

            verification(testCase, actSolution,expSolution,V)
        end
    end
end

function [] = verification(testCase, actSolution,expSolution,V)
       verifyEqual(testCase, actSolution, expSolution,'RelTol',0.2,'AbsTol',0.000000001,...
            'Difference between actual and expected exceeds relative tolerance');
       verifyLessThanOrEqual(testCase, max(V'*V), 1.02,...
            'Orthonormalization failed');
       verifyLessThanOrEqual(testCase, max(imag(V)), 0, ...
            'V is not purely real'); 
       verifyEqual(testCase, rank(V), size(V,2),...
            'Rank(V) is not full');
       verifyEqual(testCase, nnz(isinf(V)), 0, ...
            'V contains Inf');
       verifyEqual(testCase, nnz(isnan(V)), 0, ...
            'V contains Nan');
end