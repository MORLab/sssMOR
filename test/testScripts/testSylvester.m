classdef testSylvester< matlab.unittest.TestCase
    
    properties
         sysCell;
         testPath;
    end
 
    methods(TestMethodSetup)
        function getBenchmarks(testCase)
            testCase.testPath=pwd;
            if exist('benchmarksSysCell.mat','file')
                load('benchmarksSysCell.mat');
                testCase.sysCell=benchmarksSysCell;
            end
            
            %the directory "benchmark" is in sssMOR
            p = mfilename('fullpath'); k = strfind(p, 'test\'); 
            pathBenchmarks = [p(1:k-1),'benchmarks'];
            cd(pathBenchmarks);
        end
    end
    
    methods(TestMethodTeardown)
        function changePath(testCase)
            cd(testCase.testPath);
        end
    end
    
    %Test functions
    methods (Test)  
        function testSylvesterV(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                badBenchmarks = {'LF10.mat','beam.mat','random.mat'};
                if ~any(strcmp(sys.Name,badBenchmarks))
                
                % get irka shifts and tangential directions
                n = 10; r = ones(sys.m,n); l = ones(sys.p,n);
                sysrIrka = irka(sys, zeros(1,n),r,l);
                
                Opts.rType = 'dir';[r,p] = residue (sysrIrka,Opts);
                s0 = -(conj(p)); l = r{1}; r = r{2}.';         
                
                [sysr, V, ~, Bb, Rsylv] = rk(sys,s0,s0,r,l);
                warning off
                [R, S, B_] = getSylvester(sys, sysr, V);
                warning on
                
                %   compute the residuals
                res0 = norm(sysr.A - sysr.E*S - sysr.B*R);
                res1 = norm(sys.A*V - sys.E*V*S - sys.B*R);
                res2 = norm(sys.A*V - sys.E*V*(sysr.E\sysr.A)-B_*R);
                eigdiff = norm(cplxpair(eig(S))-cplxpair(s0.'));
                
%                 S2 = sysr.e\(sysr.a - sysr.b*Rsylv);
%                 eigdiff = norm(cplxpair(eig(S2))-cplxpair(s0.'))
%                 res0 = norm(sysr.A - sysr.E*S2 - sysr.B*Rsylv)
%                 res1 = norm(sys.A*V - sys.E*V*S2 - sys.B*Rsylv)
%                 res2 = norm(sys.A*V - sys.E*V*(sysr.E\sysr.A)-Bb*Rsylv)
                
                actSolution={B_, R, res0, res1, res2,eigdiff};
                expSolution={Bb, Rsylv, 0, 0, 0, 0};
                
                verification (testCase, actSolution,expSolution);
                end
            end
        end
        
        function testSylvesterW(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                badBenchmarks = {'LF10.mat','beam.mat','random.mat'};
                if ~any(strcmp(sys.Name,badBenchmarks))
                
                % get irka shifts and tangential directions
                n = 10; r = ones(sys.m,n); l = ones(sys.p,n);
                sysrIrka = irka(sys, zeros(1,n),r,l);
                
                Opts.rType = 'dir';[r,p] = residue (sysrIrka,Opts);
                s0 = -(conj(p)); l = r{1}; r = r{2}.';         
                
                [sysr, ~, W, ~, ~, Cb, Lsylv] = rk(sys,s0,s0,r,l);
                warning off
                [L, S, C_] = getSylvester(sys, sysr, W, 'W');
                warning on
                
                %   get dual system
                sys = sys.'; sysr = sysr.';
                
                %   compute the residuals
                res0 = norm(sysr.A - sysr.E*S.' - sysr.B*L);
                res1 = norm(sys.A*W - sys.E*W*S.' - sys.B*L);
                res2 = norm(sys.A*W - sys.E*W*(sysr.E\sysr.A)-C_.'*L);
                eigdiff = norm(cplxpair(eig(S))-cplxpair(s0.'));
                
%                 S2 = sysr.e\(sysr.a - sysr.b*Rsylv);
%                 eigdiff = norm(cplxpair(eig(S2))-cplxpair(s0.'))
%                 res0 = norm(sysr.A - sysr.E*S2 - sysr.B*Rsylv)
%                 res1 = norm(sys.A*V - sys.E*V*S2 - sys.B*Rsylv)
%                 res2 = norm(sys.A*V - sys.E*V*(sysr.E\sysr.A)-Bb*Rsylv)
                
                actSolution={C_, L, res0, res1, res2,eigdiff};
                expSolution={Cb, Lsylv, 0, 0, 0, 0};
                
                verification (testCase, actSolution,expSolution);
                end
            end
        end
    end
end
    
function [] = verification (testCase, actSolution,expSolution)

          % SEEMS LIKE NUMERICS IS EATING UP A LOT OF DIGITS....
          
          % do the results from rk match with those of getSylvester?
          verifyEqual(testCase, actSolution{1},expSolution{1},'RelTol',1e-6,...
              'rk and getSylvester do not match')
          verifyEqual(testCase, actSolution{2},expSolution{2},'RelTol',1e-6,...
              'rk and getSylvester do not match')
          % are the residual zero?
            verifyEqual(testCase, actSolution{3},expSolution{3},'AbsTol',1e-6,...
                'residual is not 0')          
            verifyEqual(testCase, actSolution{4},expSolution{4},'AbsTol',1e-6,...
                'residual is not 0') 
          % are the eigenvalues of S the shifts s0?
            verifyEqual(testCase, actSolution{5},expSolution{5},'AbsTol',1e-6,...
                'eigenvalues of S do not match')  
end