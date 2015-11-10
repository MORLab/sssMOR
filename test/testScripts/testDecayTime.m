classdef testDecayTime < matlab.unittest.TestCase
    
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

            for i=1:length(files)
                testCase.sysCell{i} = loadSss(files(i).name);
            end

            % change path back
            cd(Path);
        end
    end
    
    %Test functions
    methods (Test)  
        function testDecayTime1(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                
                [actTmax]=decayTime(sys);
                pmin=-log(100)/actTmax; %real(smallest dominant pole)
                [h]=impulse(ss(pmin,1,1,0),0:actTmax:actTmax); %impulse answer of pmin
                h=abs(h);
                expSolution=h(1)*0.01;
                actSolution=h(2);
                
                verification (testCase, actSolution, expSolution);
                verifyInstanceOf(testCase, actTmax , 'double', 'Instances not matching');
                verifySize(testCase, actTmax, [1 1], 'Size not matching');
            end
        end
    end
end

function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.1,'AbsTol',0.000001, ...
               'Difference between actual and expected exceeds relative tolerance');
end