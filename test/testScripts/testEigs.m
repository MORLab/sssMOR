classdef testEigs < matlab.unittest.TestCase
    
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
        function testEigs1(testCase)
            for i=1:length(testCase.sysCell)
                sys_sss=testCase.sysCell{i};
                
                [actD]=eigs(sys_sss,12);
                [expD]=eigs(full(sys_sss.A), full(sys_sss.E),12);
                
                %sort eigenvalues on real part
                actD_real=real(actD);
                actD_imag=imag(actD);
                tbl=table(actD_real, actD_imag);
                tbl=sortrows(tbl);
                actD_real=tbl.actD_real;
                actD_imag=tbl.actD_imag;

                expD_real=real(expD);
                expD_imag=imag(expD);
                tbl=table(expD_real, expD_imag);
                tbl=sortrows(tbl);
                expD_real=tbl.expD_real;
                expD_imag=tbl.expD_imag;
                
                actSolution={full(actD_real), sort(full(actD_imag))};
                expSolution={expD_real, sort(expD_imag)};
                
                verification (testCase, actSolution, expSolution);
                verifyInstanceOf(testCase, actD , 'double', 'Instances not matching');
                verifySize(testCase, actD, size(expD), 'Size not matching');
            end
        end
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.2,'AbsTol',0.00001, ...
               'Difference between actual and expected exceeds relative tolerance');
end