classdef testIssd < matlab.unittest.TestCase
    
    properties 
        pwdPath
        sysCell
        deleteBenchmarks
        testPath
    end

    methods(TestClassSetup)
        function getBenchmarks(testCase)
            testCase.pwdPath=pwd;
            if exist('benchmarksSysCell.mat','file')
                testCase.deleteBenchmarks=0;
            else
                testCase.testPath=loadBenchmarks;
                testCase.deleteBenchmarks=1;
            end
            
            temp=load('benchmarksSysCell.mat');
            testCase.sysCell=temp.benchmarksSysCell;
            if isempty(testCase.sysCell)
                error('No benchmarks loaded.');
            end

            %the directory "benchmark" is in sssMOR
            p = mfilename('fullpath'); k = strfind(p, 'test\'); 
            pathBenchmarks = [p(1:k-1),'benchmarks'];
            cd(pathBenchmarks);
        end
    end
    
    methods(TestClassTeardown)
        function changePath(testCase)
            if testCase.deleteBenchmarks
                cd(testCase.testPath);
                delete('benchmarksSysCell.mat');
            end
            cd(testCase.pwdPath);
        end
    end
    
    %Test functions
    methods (Test)  
        function testIssd1(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                
                size(sys.A)
                [actSolution]=issd(sys);
                if nnz(eig(full(sys.A+sys.A'),'vector')>=0)
                    expSolution=0;
                else
                    expSolution=1;
                end
                if sys.isDescriptor
                    if nnz(eig(full(sys.E))<=0)
                        expSolution=0;
                    end
                end

                verification (testCase, actSolution, expSolution);
                verifyInstanceOf(testCase, actSolution , 'double', 'Instances not matching');
                verifySize(testCase, actSolution, [1 1], 'Size not matching');
            end
        end
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.1,'AbsTol',0.000001, ...
               'Difference between actual and expected exceeds relative tolerance');
end