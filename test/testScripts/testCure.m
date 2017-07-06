classdef testCure < sssTest
   
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
                warning off, sysr = cure(sys); warning on 
                
                verification(testCase, {}, {}, sysr, sys);
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
                          'stopval',5e-2,...
                          'verbose', 1,...
                          'test', 1,...
                          'gif',1,...
                          'maxIter',10);
            Opts.warn =1;
            Opts.w = w;

            warning off, sysr = cure(sys,Opts); warning on
            close all, delete('CURE.gif')
            
            verification(testCase, {}, {}, sysr, sys);
            verifyLessThanOrEqual(testCase, sysr.n,Opts.cure.maxIter*2);
        end
        function testCuredIrka(testCase)
            sys = loadSss('rail_1357');
            disp(sys)
            if ~sys.isSiso, sys = sys(1,1); end

            % run cured spark with non default values
            Opts.cure = struct('redfun','irka',...
                          'fact','V',...
                          'init','sm',...
                          'stop','nmax',...
                          'stopval',8,...
                          'nk', 4);
            warning off
            sysr = cure(sys,Opts);
            verifyEqual(testCase, sysr.n, Opts.cure.stopval,'Stopcrit violated');

            
            Opts.cure = struct('redfun','irka',...
                          'fact','W',...
                          'init','slm',...
                          'stop','nmax',...
                          'stopval',8,...
                          'nk', 4);
            sysr = cure(sys,Opts);
            warning on
            
            verifyEqual(testCase, sysr.n, Opts.cure.stopval,'Stopcrit violated');
        end
        function testCuredRk(testCase)
            sys = loadSss('eady');
            disp(sys)
            if ~sys.isSiso, sys = sys(1,1); end
            warning off
            % run cured spark with non default values
            Opts.cure = struct('redfun','rk+pork',...
                          'fact','V',...
                          'init','slm');
            cure(sys,Opts);
            Opts.cure = struct('redfun','rk+pork',...
                          'fact','W',...
                          'init','slm');
            sysr = cure(sys,Opts);
            warning on
            
            verification(testCase, {}, {}, sysr, sys);
        end
    end
end

function [] = verification(testCase, ~, ~, sysr, sys)
       verifyLessThanOrEqual(testCase, sysr.n, sys.n,'ROM bigger than ROM'); 
       verifyTrue(testCase,isstable(sysr),'ROM unstable');
end