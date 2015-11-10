classdef testEig < matlab.unittest.TestCase
    
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
            warning off
            for i=1:length(files)
                testCase.sysCell{i} = loadSss(files(i).name);
            end
            warning on

            % change path back
            cd(Path);
        end
    end
    
    %Test functions
    methods (Test)  
        function testEig1(testCase)
            for i=1:length(testCase.sysCell)
                sys_sss=testCase.sysCell{i};
                
                [actD]=eig(sys_sss);
                [expD]=eig(full(sys_sss.A), full(sys_sss.E));
                
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
                
                actSolution={full(actD_real), full(actD_imag)};
                expSolution={expD_real, expD_imag};
                
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