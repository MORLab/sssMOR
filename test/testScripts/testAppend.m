classdef testAppend < matlab.unittest.TestCase
    
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
            function testAppend1(testCase)
            for i=1:2:length(testCase.sysCell)-1
                sys1_sss=testCase.sysCell{i};
                if sys1_sss.isDescriptor
                    sys1_ss=ss(full(sys1_sss.E\sys1_sss.A),full(sys1_sss.E\sys1_sss.B),full(sys1_sss.C),0);
                else
                    sys1_ss=ss(sys1_sss);
                end
                sys2_sss=testCase.sysCell{i+1};
                if sys2_sss.isDescriptor
                    sys2_ss=ss(full(sys2_sss.E\sys2_sss.A),full(sys2_sss.E\sys2_sss.B),full(sys2_sss.C),0);
                else
                    sys2_ss=ss(sys2_sss);
                end

                expSys = append(sys1_ss, sys2_ss);
                actSys = append(sys1_sss, sys2_sss);

                expSolution={expSys.A, expSys.B, expSys.C, expSys.D};
                actSolution={full(actSys.E\actSys.A), full(actSys.E\actSys.B), full(actSys.C), full(actSys.D)};

                verification (testCase, actSolution, expSolution);
                verifyInstanceOf(testCase, actSys.A , 'double', 'Instances not matching');
                verifyInstanceOf(testCase, actSys.B , 'double', 'Instances not matching');
                verifyInstanceOf(testCase, actSys.C , 'double', 'Instances not matching');
                verifyInstanceOf(testCase, actSys.D , 'double', 'Instances not matching');
                verifyInstanceOf(testCase, actSys.E , 'double', 'Instances not matching');
                verifySize(testCase, full(actSys.A), size(expSys.A), 'Size not matching');
                verifySize(testCase, full(actSys.B), size(expSys.B), 'Size not matching');
                verifySize(testCase, full(actSys.C), size(expSys.C), 'Size not matching');
                verifySize(testCase, full(actSys.D), size(expSys.D), 'Size not matching');
                verifySize(testCase, full(actSys.E), size(actSys.A), 'Size not matching');
            end
        end 
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.1,'AbsTol',0.000001, ...
               'Difference between actual and expected exceeds relative tolerance');
end