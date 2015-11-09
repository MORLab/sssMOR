classdef testPork < matlab.unittest.TestCase
    
    properties
         sysCell;
    end
 
    methods(TestMethodSetup)
        function getBenchmarks(testCase)
            % change path
            Path = pwd; %original
            
            %insert path of local benchmark folder
            %the directory "benchmark" is in sssMOR
            p = mfilename('fullpath'); k = strfind(p, 'test\'); 
            pathBenchmarks = [p(1:k-1),'benchmarks'];
            cd(pathBenchmarks);

            % load files
            files = dir('*.mat'); 
            testCase.sysCell=cell(1,length(files));
%             warning('off','sssMOR:loadSss:2ndOrder')
            for i=1:length(files)
                testCase.sysCell{i} = loadSss(files(i).name);
            end
%             warning('on','loadSss:2ndOrder')

            % change path back
            cd(Path);
        end
    end
    
    %Test functions
    methods (Test)  
        function testPorkV(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                
                n = 10; s0 = rand(1,n);
                [sysr, V] = rk(sys,s0);
                [Crt, S_V] = getSylvester(sys, sysr, V);
                
                [Ar,Br,Cr,Er] = porkV(V,S_V,Crt,sys.C);
                sysr = sss(Ar,Br,Cr,sys.D,Er);
                
                actSolution={sysr,s0,sys};
                
                verification (testCase, actSolution);
            end
        end
    end
end
    
function [] = verification (testCase, actSolution)
          sysr = actSolution{1}; s0 = actSolution{2}; sys = actSolution{3};
          % are all matrices double
          verifyInstanceOf(testCase, sysr.a , 'double', 'Instances not matching');
          verifyInstanceOf(testCase, sysr.b , 'double', 'Instances not matching');
          verifyInstanceOf(testCase, sysr.c , 'double', 'Instances not matching');
          verifyInstanceOf(testCase, sysr.d , 'double', 'Instances not matching');
          verifyInstanceOf(testCase, sysr.e , 'double', 'Instances not matching');
          % are all matrices of the correct size?
          verifySize(testCase, sysr.a, [length(s0),length(s0)], 'Size not matching');
          % is Er == I?
          verifyEqual(testCase, sysr.e , eye(length(s0)), 'Er~= I');
          % does Ar have eigenvalues at the mirror images of the shifts?
           verifyEqual(testCase, eig(sysr.a) , -s0 , 'Eigenvalues do not match');
           % do the moments match?
           for iS = 1:length(s0)
           verifyEqual(testCase, moments(sys,s0(iS),1) , moments(sysr,s0(iS),1), 'moments do not match');
           end
end