classdef testSolveLse < sssTest
    methods(Test)
        function testImpulseBasic(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if sys.isSiso && ~sys.isDae
                    % test A\B
                    actX1=solveLse(sys.A,sys.B);
                    expX1=full(sys.A\sys.B);
                    
                    % test A\B and A'\C'
                    Opts.lse='hess';
                    [actX2,actY2]=solveLse(sys.A,sys.B,sys.C,Opts);
                    expX2=full(sys.A\sys.B);
                    expY2=full(sys.A'\sys.C');
                    
                    % test (A-5*E)\B and (A-5*E)'\C'
                    Opts.lse='hess';
                    [actX3,actY3]=solveLse(sys.A,sys.B,sys.C,sys.E,[5,8,3],Opts);
                    expX3=[full((sys.A-5*sys.E)\sys.B),full((sys.A-8*sys.E)\sys.B),full((sys.A-3*sys.E)\sys.B)];
                    expY3=[full((sys.A-5*sys.E)'\sys.C'),full((sys.A-8*sys.E)'\sys.C'),full((sys.A-3*sys.E)'\sys.C')];
                    
                    % test Markov parameter E\B
                    [actX4]=solveLse(sys.A,sys.B,sys.E,Inf);
                    expX4=full(sys.E\sys.B);
                    
                    % test several shifts
                    [actX5]=solveLse(sys.A,sys.B,sys.E,[3,4,4, 1+i, 1-i]);
                    expX5=full([(sys.A-3*sys.E)\sys.B,(sys.A-4*sys.E)\sys.B,...
                        (sys.A-4*sys.E)\sys.B,(sys.A-(1+i)*sys.E)\sys.B,...
                        (sys.A-(1-i)*sys.E)\sys.B]);
                    
                    % test iterative solver
                    Opts.lse='iterative';
                    actX6=solveLse(sys.A,sys.B,sys.E,6,Opts);
                    expX6=full((sys.A-6*sys.E)\sys.B);
                    
                    actSolution={actX1,actX2,actY2,actX3,actY3,actX4,actX5,actX6};
                    expSolution={expX1,expX2,expY2,expX3,expY3,expX4,expX5,expX6};
                    
                	verification(testCase, actSolution, expSolution);
                end
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1e-5,'AbsTol',1e-5,...
    'Difference between actual and expected exceeds relative tolerance');
end