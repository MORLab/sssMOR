classdef testImpulse < matlab.unittest.TestCase
    
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
        function testImpulse1(testCase)
            for i=1:length(testCase.sysCell)
                sys_sss=testCase.sysCell{i};
                sys_ss=ss(sys_sss);
                
                [actH, actT]=impulse(sys_sss); 
                [expH]=impulse(ss(sys_ss), actT);
                         
                
                actSolution={actH{1}'};
                expSolution={expH};
                
                verification (testCase, actSolution, expSolution);
                %returns cell
%                 verifyInstanceOf(testCase, actT , 'double', 'Instances not matching');
%                 verifyInstanceOf(testCase, actH , 'double', 'Instances not matching');
%                 verifySize(testCase, actH, size(expH), 'Size not matching');
%                 verifySize(testCase, actT, size(expT), 'Size not matching');
            end
        end
%         function testImpulse2(testCase)
%             %does not work
%             for i=1:length(testCase.sysCell)
%                 sys_sss=testCase.sysCell{i};
%                 sys_ss=ss(sys_sss);
%                 t=1:0.1:5;
%                 
%                 [actH]=impulse(sys_sss,t); 
%                 [expH]=impulse(ss(sys_ss), t);
%                          
%                 
%                 actSolution={actH{1}'};
%                 expSolution={expH};
%                 
%                 verification (testCase, actSolution, expSolution);
% %                 verifyInstanceOf(testCase, actT , 'double', 'Instances not matching');
% %                 verifyInstanceOf(testCase, actH , 'double', 'Instances not matching');
% %                 verifySize(testCase, actH, size(expH), 'Size not matching');
% %                 verifySize(testCase, actT, size(expT), 'Size not matching');
%             end
%         end
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.1,'AbsTol',0.000001, ...
               'Difference between actual and expected exceeds relative tolerance');
end