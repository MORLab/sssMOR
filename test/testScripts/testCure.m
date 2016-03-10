classdef testCure < matlab.unittest.TestCase
    
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
        function testCureDef(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                badBenchmarks = {'LF10.mat','beam.mat','random.mat',...
                    'SpiralInductorPeec.mat'};
                
                if ~any(strcmp(sys.Name,badBenchmarks))
                disp(sys)
                if ~sys.isSiso, sys = sys(1,1); end
                
                % run cure (default options)
                warning off
                cure(sys);
                warning on 
                
                % no verification available for this
                end
            end
        end
        function testCureNonDef(testCase)
            sys = loadSss('iss');
            disp(sys)
            if ~sys.isSiso, sys = sys(1,1); end
            [~, w] = freqresp(sys);

            % run cured spark with non default values
            Opts.cure = struct('fact','W',...
                          'init','slm',...
                          'stop','h2Error',...
                          'stopval',norm(sys)*1e-1,...
                          'verbose', 1,...
                          'test', 1,...
                          'gif',1,...
                          'maxIter',10);
            Opts.warn =1;
            Opts.w = w;

            warning off
            cure(sys,Opts);
            warning on
            close all
        end
        function testCuredIrka(testCase)
            sys = loadSss('rail_1357');
            disp(sys)
            if ~sys.isSiso, sys = sys(1,1); end

            % run cured spark with non default values
            Opts.cure = struct('redfun','irka',...
                          'fact','V',...
                          'init','slm',...
                          'stopval',20,...
                          'nk', 4);
            warning off
            cure(sys,Opts);
            
            Opts.cure = struct('redfun','irka',...
                          'fact','W',...
                          'init','slm',...
                          'stopval',20,...
                          'nk', 4);
            cure(sys,Opts);
            warning on
        end
        function testCuredRk(testCase)
            sys = loadSss('eady');
            disp(sys)
            if ~sys.isSiso, sys = sys(1,1); end
            warning off
            % run cured spark with non default values
            Opts.cure = struct('redfun','rk+pork',...
                          'fact','V',...
                          'init','slm',...
                          'stopval',20,...
                          'nk', 5);
            cure(sys,Opts);
            Opts.cure = struct('redfun','rk+pork',...
                          'fact','W',...
                          'init','slm',...
                          'stopval',20,...
                          'nk', 5);
            cure(sys,Opts);
            warning on
        end
    end
end