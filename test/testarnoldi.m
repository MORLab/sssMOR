classdef testarnoldi < matlab.unittest.TestCase

 
    properties
    end
 
    methods(Test)     
        
        function testarnoldi1(testCase) 
            load build;

            [V,AV,EV] = arnoldi(speye(size(A)),A,B,ones(1,10)*100);
            actSolution={V(:,1),AV(:,1),EV(:,1)};

            temp=(A-100*eye(size(A)))\B;
            expV=temp/norm(temp);
            expSolution={expV,A*expV,eye(size(A))*expV};

            verification(testCase, actSolution,expSolution,V)
        end
        
        function testarnoldi2(testCase) 
            load beam;

            [V,AV,EV] = arnoldi(speye(size(A)),A,B,[ones(1,5)*100*1i,-ones(1,5)*100*1i]);
            actSolution={V(:,1),AV(:,1),EV(:,1)};

            temp=real((A-100*1i*eye(size(A)))\B);
            expV=temp/norm(temp);
            expSolution={expV,A*expV,eye(size(A))*expV};

            verification(testCase, actSolution,expSolution,V)

        end
        
        function testarnoldi3(testCase)
            load fom;

            [V,AV,EV] = arnoldi(speye(size(A)),A,B,1:10);
            actSolution={V(:,1),AV(:,1),EV(:,1)};

            temp=(A-1*eye(size(A)))\B;
            expV=temp/norm(temp);
            expSolution={expV,A*expV,eye(size(A))*expV};

            verification(testCase, actSolution,expSolution,V)
            
        end
        
        function testarnoldi4(testCase) 
            load eady;

            [V,AV,EV] = arnoldi(speye(size(A)),A,B,[1+1i, 1-1i, 2+2i, 2-2i, 3+3i, 3-3i, 4+4i, 4-4i, 5+5i, 5-5i]);
            actSolution={V(:,1),AV(:,1),EV(:,1)};

            temp=real((A-(1+1i)*eye(size(A)))\B);
            expV=temp/norm(temp);
            expSolution={expV,A*expV,eye(size(A))*expV};
            
            verification(testCase, actSolution,expSolution,V)

        end
        
        function testarnoldi5(testCase) 
            load random;
            [V,AV,EV] = arnoldi(speye(size(A)),A,B, [Inf,  Inf]);
            actSolution={V(:,1),AV(:,1),EV(:,1)};

            temp=full(real(B));
            expV=temp/norm(temp);
            expSolution={expV,A*expV,eye(size(A))*expV};
            
            verifyEqual(testCase, actSolution, expSolution,'RelTol',0.1,...
            'Difference between actual and expected exceeds relative tolerance');

        end
        
    end
end

function [] = verification(testCase, actSolution,expSolution,V)
       verifyEqual(testCase, actSolution, expSolution,'RelTol',0.1,...
            'Difference between actual and expected exceeds relative tolerance');
       verifyLessThanOrEqual(testCase, max(V'*V), 1.02,...
            'Orthonormalization failed');
       verifyLessThanOrEqual(testCase, max(imag(V)), 0, ...
            'V is not purely real'); 
       verifyEqual(testCase, rank(V), 10,...
            'Rank(V) is not full');
       verifyEqual(testCase, nnz(isinf(V)), 0, ...
            'V contains Inf');
       verifyEqual(testCase, nnz(isnan(V)), 0, ...
            'V contains Nan');
end