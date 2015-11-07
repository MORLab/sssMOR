classdef testBode < matlab.unittest.TestCase
    
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
        function testBode1(testCase)
            for i=1:length(testCase.sysCell)
                sys_sss=testCase.sysCell{i};
                sys_ss=ss(sys_sss);
                t=1:100:1000;
                
                [actMag, actPhase, actOmega]=bode(sys_sss,t);
                [expMag, expPhase, expOmega]=bode(ss(sys_ss),t);
                    
                %Phase between 0° to 360°
                for j=1:length(actPhase)
                    if actPhase(:,:,j)<0
                    actPhase(:,:,j)=actPhase(:,:,j)+360;
                    end
                    if expPhase(:,:,j)<0
                    expPhase(:,:,j)=expPhase(:,:,j)+360;
                    end
                end
                
                actSolution={actMag, actPhase, actOmega};
                expSolution={expMag, expPhase, expOmega};
                
                verification (testCase, actSolution, expSolution);
                verifyInstanceOf(testCase, actMag , 'double', 'Instances not matching');
                verifyInstanceOf(testCase, actPhase , 'double', 'Instances not matching');
                verifyInstanceOf(testCase, actOmega , 'double', 'Instances not matching');
                verifySize(testCase, actMag, size(expMag), 'Size not matching');
                verifySize(testCase, actPhase, size(expPhase), 'Size not matching');
                verifySize(testCase, actOmega, size(expOmega), 'Size not matching');
            end
        end
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.1,'AbsTol',0.000001, ...
               'Difference between actual and expected exceeds relative tolerance');
end