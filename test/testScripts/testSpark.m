classdef testSpark < sssTest
    
    methods (Test)  
        function testSparkDef(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                badBenchmarks = {'LF10.mat','beam.mat','random.mat'};
                
                if ~any(strcmp(sys.Name,badBenchmarks))
                disp(sys)
                if ~sys.isSiso, sys = sys(1,1); end
                
                % run spark (default options)
                s0 = rand(1,2);
                [V,Sv,Rv] = spark(sys,s0);
                
                % check Sylvester equation
                res = norm(sys.A*V - sys.E*V*Sv - sys.B*Rv);
                
                actSolution={V,Sv,Rv,res,sys};
                
                verification (testCase, actSolution);
                end
            end
        end
        function testSparkNonDef(testCase)
            sys = loadSss('iss');
            disp(sys)
            if ~sys.isSiso, sys = sys(1,1); end

            % run spark (default options)
            s0 = rand(1,2);
            Opts.spark = struct('type','standard',...
                          'test',1,...
                          'verbose',1,...
                          'mfe', 5e2,...
                          'mi', 50,...
                          'xTol',1e-5,...
                          'fTol',1e-5,...
                          'modelTol',1e-3);

            [V,Sv,Rv] = spark(sys,s0,Opts);
            close all
            % check Sylvester equation
            res = norm(sys.A*V - sys.E*V*Sv - sys.B*Rv);

            actSolution={V,Sv,Rv,res,sys};

            verification (testCase, actSolution);
        end
        function testSparkSysr(testCase)
            sys = loadSss('beam.mat');
            disp(sys)
            
            s0 = rand(1,2);
            
            %Reduction with porkV           
            Opts.spark.pork = 'V';
            [sysr,V,Sv,Rv,~] = spark(sys,s0,Opts);
            sysrPO = porkV(V,Sv,Rv,sys.c);
            
            verifyEqual(testCase,sysr.a,sysrPO.a,'AbsTol', 1e-8, ...
                'Verification with porkV failed');
            verifyEqual(testCase,sysr.b,sysrPO.b,'AbsTol', 1e-8, ...
                'Verification with porkV failed');
            verifyEqual(testCase,sysr.c,sysrPO.c,'AbsTol', 1e-8, ...
                'Verification with porkV failed');
            verifyEqual(testCase,sysr.e,sysrPO.e,'AbsTol', 1e-8, ...
                'Verification with porkV failed');
            
            %Reduction with porkW
            Opts.spark.pork = 'W';
            [sysr,W,Sw,Rw,~] = spark(sys.',s0,Opts);
            Sw = Sw';
            sysrPO = porkW(W,Sw,Rw,sys.b);
            
            verifyEqual(testCase,sysr.a,sysrPO.a,'AbsTol', 1e-8, ...
                'Verification with porkW failed');
            verifyEqual(testCase,sysr.b,sysrPO.b,'AbsTol', 1e-8, ...
                'Verification with porkW failed');
            verifyEqual(testCase,sysr.c,sysrPO.c,'AbsTol', 1e-8, ...
                'Verification with porkW failed');
            verifyEqual(testCase,sysr.e,sysrPO.e,'AbsTol', 1e-8, ...
                'Verification with porkW failed');
        end
    end
end
    
function [] = verification (testCase, actSolution)
          V = actSolution{1}; S = actSolution{2}; 
          R = actSolution{3}; res = actSolution{4}; sys = actSolution{5};
          % sizes
          verifySize(testCase, V, [sys.n, 2], 'Size of V not matching');
          verifySize(testCase, S, [2,2] , 'Size of S not matching');
          verifySize(testCase, R, [sys.m,2] , 'Size of R not matching');
          % Check residual of Sylvester equation
          verifyEqual(testCase, res, 0, 'AbsTol', 1e-8, ...
              'Sylvester equation not satisfied');
end