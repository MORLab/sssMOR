classdef testAllOpts < sssTest
% testAllOpts - testing of all approaches with all possible combination of
%               parameters
%
% Description:
%  
% Authors:      
% Last Change:  25 Oct 2017
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------
 
    methods(Test)
        % Extras
        function runInitializeShifts(testCase)
            % Define possible opts combinations
            OptsStuct.strategy   = {'eigs','ADI','const', 'ROM', 'linspaced', 'logspaced', 'random', 'lognrnd'};
            OptsStuct.shiftType  = {'conj','real', 'imag'};
            OptsStuct.wmin       = {abs(eigs(sys,1,'sm'))}; %TODO any way of predefining the wmin, wmax params without sys?
            OptsStuct.wmax       = {abs(eigs(sys,1))};
            OptsStuct.kp         = {40};
            OptsStuct.km         = {25};
            OptsStuct.eigsType   = {'sm', 'lm', 'li', 'si', 'lr', 'sr', 'la','sa'};
            OptsStuct.constValue = {0};
            OptsStuct.offset     = {0};
            OptsStuct.format     = {'complex','ab'};

            [AllOptsCell,nCases]  = generateAllOpts(OptsStuct);

            h = waitbar(0,'initializeShifts: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test
                    sys = testCase.sysCell{i};
                    initializeShifts(sys,Opts);
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)                    
            end
        end
        function runIsH2opt(testCase)
            % Define possible opts combinations
            OptsStuct.tol   = {1e-3};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStuct);

            h = waitbar(0,'isH2opt: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test
                    sys = testCase.sysCell{i};
                    [sysr, ~, ~, s0opt]  = irka(sys, -eigs(sys,8).');
                    isH2opt(sys, sysr, s0opt);
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)                    
            end
        end
        function runMoments(testCase)
            % Define possible opts combinations
            OptsStuct.refine      = {0};
            OptsStuct.refTol      = {1e-6};
            OptsStuct.refMaxiter  = {1e2};
            OptsStuct.lse         = {'sparse'};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStuct);

            h = waitbar(0,'moments: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test
                    sys = testCase.sysCell{i};
                    s0=5;
                    moments(sys, s0, 4, Opts);
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)                    
            end
        end

        % Classic
        function runArnoldi(testCase)
            % Define possible opts combinations
            OptsStuct.real        = {true, false};
            OptsStuct.orth        = {false, 'dgks', 'mgs'};
            OptsStuct.reorth      = {'mgs', false, 'qr'};
            OptsStuct.lse         = {'sparse','full','hess','iterative'};
            OptsStuct.dgksTol     = {1e-12};
            OptsStuct.krylov      = {'standardKrylov','cascadedKrylov'};
            OptsStuct.dgksTol     = {1e3};
            OptsStuct.tollse      = {1e-6};
            OptsStuct.solver      = {'cgs','bicgstab','bicg'};
            OptsStuct.verbose     = {true, false};
            OptsStuct.force       = {true, false};

            [AllOptsCell,nCases]  = generateAllOpts(OptsStuct);

        %     verifyError(testCase,arnoldi(),'sss:solveLse:cascadeSiso')

            h = waitbar(0,'arnodi: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test
                    sys = testCase.sysCell{i};
                    sysr = arnoldi(speye(size(sys.A)),sys.A,sys.B,[1-1i, 1-1i, 1-1i, 1+1i, 1+1i, 1+1i], Opts);
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)                    
            end
        end
        function runModalMor(testCase)
            % Define possible opts combinations
            OptsStuct.type        = {'SM','LM',1};
            OptsStuct.orth        = {true, false,'qr'};
            OptsStuct.real        = {true, false,'real'};
            OptsStuct.tol         = {1e-6};
            OptsStuct.dominance   = {0,'analyze','2q','3q'};
            OptsStuct.lse         = {'sparse','full'};
            OptsStuct.subspaceW   = {'eigs','1by1'};

            [AllOptsCell,nCases]  = generateAllOpts(OptsStuct);

            h = waitbar(0,'modalMor: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test
                    sys = testCase.sysCell{i};
                    if strcmp(Opts.dominance,'2q') && strcmp(testCase.sysCell{i}.Name,'iss')
                        % Option dominance = '2q' produces a numerical error
                        % for the benchmark "iss", because only 2 of the 
                        % requested 8 eigenvalues converge.

                        warning('off')
                        try
                            sysr = modalMor(sys,4, Opts);
                        end
                        warning('on')
                    else
                        sysr = modalMor(sys,4, Opts);
                    end
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)                    
            end
        end
        function runProjectiveMor(testCase)
            % Define possible opts combinations
            OptsStuct.trans       = {'T','H'};

            [AllOptsCell,nCases]  = generateAllOpts(OptsStuct);

            h = waitbar(0,'projectiveMor: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test

                    sys = testCase.sysCell{1};
                    s0 = 1; r = ones(sys.m,length(s0));

                    V = (sys.A-s0*sys.E)\(sys.B*r);
                    sysr = projectiveMor(sys,V,Opts);
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)                    
            end
        end
        function runRk(testCase)
            % Define possible opts combinations
            OptsStuct.real       = {true,false};

            [AllOptsCell,nCases]  = generateAllOpts(OptsStuct);

            h = waitbar(0,'rk: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test
                    sys = testCase.sysCell{1};
                    n = 5; s0val = 100; s0 = [ones(1,n)*s0val*1i,-ones(1,n)*s0val*1i]; 
                    [sysr, ~] = rk(sys, [], s0, Opts);
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)
            end
        end
        function runTbr(testCase)
            % Define possible opts combinations
            OptsStuct.type        = {'tbr', 'adi', 'matchDcGain'};
            OptsStuct.redErr      = {'0', 1};
            OptsStuct.hsvTol      = {1e-15};
            OptsStuct.warnOrError = {'warn','error', 0};
            OptsStuct.lse         = {'sparse', 'full','iterative', 'gauss'}; 
            OptsStuct.rctol       = {1e-15};
            OptsStuct.forceOrder  = {true, false};


            [AllOptsCell,nCases]  = generateAllOpts(OptsStuct);

            h = waitbar(0,'tbr: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test
                    sys = testCase.sysCell{1};
                    q   = 2;
                    [sysr, ~] = tbr(sys,q, Opts);
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)
            end
        end

        % State Of The Art
        function runCirka(testCase)
            % Presets
            n = 2;
            s0 = zeros(1,n);
            % Define possible opts combinations
            OptsStuct.qm0     = {2*length(s0)};
            OptsStuct.s0m     = {[s0,s0]};
            OptsStuct.maxiter = {15};
            OptsStuct.tol     = {1e-3};
            OptsStuct.stopCrit= {'s0', 'sysr', 'sysm', 'combAny', 'combAll'};

            OptsStuct.verbose         = {true, false};
            OptsStuct.plot            = {true, false};
            OptsStuct.suppressWarn    = {true, false};
            OptsStuct.updateModel     = {'new', 'all'};
            OptsStuct.modelTol        = {1e-2};
            OptsStuct.clearInit       = {true, false};
            OptsStuct.stableModelFct  = {true, false};

            OptsStuct.irka.suppressverbose = {true, false};
            OptsStuct.irka.stopCrit        = {'combAny', 's0', 'sysr', 'combAll'};
            OptsStuct.lse                  = {'sparse', 'full', 'hess'};
            OptsStuct.irka.tol             = {1e-6};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStuct);

            h = waitbar(0,'circa: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test
                    sys = testCase.sysCell{i};
                    sys = sys(1,1);
                    sysr = cirka(sys, s0, Opts);
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)                    
            end
        end
        function runCure(testCase)

%           -.cure.maxIter: [{'20'} / positive integer]

            % Define possible opts combinations
            OptsStuct.warn            = {true, false};
            OptsStuct.w               = {[]}; % [{''} / '{wmin,wmax}' / vector of frequencies]
            OptsStuct.zeroThres       = {1e-4};
            OptsStuct.cure.verbose    = {true, false};
            OptsStuct.cure.redfun     = {'spark', 'irka', 'rk+pork'};
            OptsStuct.cure.fact       = {'V', 'W'};
            OptsStuct.cure.nk         = {2};
            OptsStuct.cure.stop       = {'normROM', 'nmax', 'h2Error'};
            OptsStuct.cure.stopval    = {1e-6};
            OptsStuct.cure.initMode   = {'zero', 'sm', 'lm', 'slm'};
            OptsStuct.cure.initN      = {5*2};
            OptsStuct.cure.fact       = {'V'};
            
            OptsStuct.cure.test       = {true, false};
            OptsStuct.cure.gif        = {true, false};
            
            OptsStuct.cure.maxIter    = {20};
            OptsStuct.cure.checkEVB   = {true};
            OptsStuct.cure.sEVBTol    = {1e-16};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStuct);

            h = waitbar(0,'cure: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test
                    sys = testCase.sysCell{i};
                    sysr = cure(sys,Opts);
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)                    
            end
        end
        function runIrka(testCase)
            % Define possible opts combinations
            OptsStuct.refine      = {0};
            OptsStuct.refTol      = {1e-6};
            OptsStuct.refMaxiter  = {1e2};
            OptsStuct.lse         = {'sparse'};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStuct);

            h = waitbar(0,'irka: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test
                    sys = testCase.sysCell{i};

                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)                    
            end
        end
        function runIsrk(testCase)
            % Define possible opts combinations
            OptsStuct.maxiter     = {50}; 
            OptsStuct.tol         = {1e-3}; 
            OptsStuct.verbose     = {true, false};
            OptsStuct.stopCrit    = {'combAny', 's0', 'sysr', 'combAll'};
            OptsStuct.suppressverbose = {true, false};
            OptsStuct.lyapchol    = {'', 'adi', 'builtIn'};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStuct);

            h = waitbar(0,'isrk: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test
                    sys = testCase.sysCell{i};
                    s0=[0,0,100,1+5i,1-5i,14-0.2i,14+0.2i, Inf, Inf];
                    sysr= isrk(sys, s0, Opts);
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)                    
            end
        end
        function runModelFct(testCase)
            % Define possible opts combinations
            OptsStuct.updateModel = {'new', 'all', 'new', 'lean'};
            OptsStuct.modelTol    = {1e-2};
            OptsStuct.plot        = {true, false};
            OptsStuct.tol         = {1e-2};

            [AllOptsCell,nCases]  = generateAllOpts(OptsStuct);

            h = waitbar(0,'modelFct: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test
                    sys = testCase.sysCell{i};
%       [sysm, s0mTot, V, W] = MODELFCT(sys,s0m,s0mTot,V,W,Opts)
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)                    
            end
        end
        function runModelFctMor(testCase)
            % Preset
            n = 2;
            s0 = zeros(1,n);
            
            % Define possible opts combinations
            OptsStuct.qm0           = {2*length(s0)};
            OptsStuct.s0m           = {[s0,s0]};
            OptsStuct.maxiter       = {15};
            OptsStuct.tol           = {1e-3};
            OptsStuct.stopcrit      = {'combAny', 's0', 'sysr', 'sysm', 'combAny', 'combAll'};
            OptsStuct.verbose       = {true, false};
            OptsStuct.plot          = {true, false};
            OptsStuct.suppressWarn  = {true, false};
            OptsStuct.updateModel   = {'new', 'all'};
            OptsStuct.modelTol      = {1e-2};
            OptsStuct.clearInit     = {true, false};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStuct);

            h = waitbar(0,'modelFctMor: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test
                    sys = testCase.sysCell{i};
                    redFct   = @(sys,s) irka(sys,s);
                    redFctOut= @(sys,s) getDesiredOutput(redFct,[1,4],sys,s);
                    sysr = modelFctMor(sys,redFctOut,s0,Opts);  
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)                    
            end
        end
        function runRkIcop(testCase)
            % Define possible opts combinations
            OptsStuct.rk      = {'twoSided', 'input', 'output'}
            OptsStuct.tol     = {1e-2};
            OptsStuct.maxIter = {100};
            OptsStuct.lse     = {'sparse', 'full', 'hess', 'gauss', 'iterative'};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStuct);

            h = waitbar(0,'rkIcop: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test
                    sys = testCase.sysCell{i};
                    s0 = 1;
                    %s0 = rkOp(sys);
                    q = 20;
                    sysr = rkIcop(sys,s0,q,Opts);
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)                    
            end
        end
        function runRkOp(testCase)
            % Define possible opts combinations
            OptsStuct.rk = {'twoSided', 'input', 'output'};
            OptsStuct.lse = {'sparse', 'full', 'hess', 'gauss', 'iterative'};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStuct);

            h = waitbar(0,'rkOp: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test
                    sys = testCase.sysCell{i};
                    q = 20;
                    sysr = rkOp(sys, q, Opts);
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)                    
            end
        end
        function runSpark(testCase)
            % Define possible opts combinations
            OptsStuct.zeroThres  = {1e-6}; 

            OptsStuct.spark.type        = {'model', 'standard'};
            OptsStuct.spark.test        = {true, false};
            OptsStuct.spark.verbose     = {true, false};
            OptsStuct.spark.mfe         = {5e3};
            OptsStuct.spark.mi          = {5e3};
            OptsStuct.spark.xTol        = {1e-10};
            OptsStuct.spark.fTol        = {1e-10};
            OptsStuct.spark.modelTol    = {1e-5};
            OptsStuct.spark.pork        = {'V', 'W'};

            OptsStuct.mespark.ritz      = {true, false};
            OptsStuct.mespark.pertIter  = {5};
            OptsStuct.mespark.maxIter   = {20};

            [AllOptsCell,nCases]        = generateAllOpts(OptsStuct);

            h = waitbar(0,'spark: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    % Run test
                    sys = testCase.sysCell{i};
                    s0 = rand(1,2);
                    [V,Sv,Rv,k] = spark(sys,s0,Opts);
                end
            end
            close(h)
            catch err
                close(h)
                sys, Opts
                fprintf(2,'Following error occurred with the options above:\n');
                rethrow(err)                    
            end
        end
    end
end