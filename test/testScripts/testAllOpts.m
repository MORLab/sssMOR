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
            OptsStruct.strategy   = {'eigs','ADI','const', 'ROM', 'linspaced', 'logspaced', 'random', 'lognrnd'};
            OptsStruct.shiftType  = {'conj','real', 'imag'};
            OptsStruct.wmin       = {abs(eigs(sys,1,'sm'))}; %TODO any way of predefining the wmin, wmax params without sys?
            OptsStruct.wmax       = {abs(eigs(sys,1))};
            OptsStruct.kp         = {40};
            OptsStruct.km         = {25};
            OptsStruct.eigsType   = {'sm', 'lm', 'li', 'si', 'lr', 'sr', 'la','sa'};
            OptsStruct.constValue = {0};
            OptsStruct.offset     = {0};
            OptsStruct.format     = {'complex','ab'};

            [AllOptsCell,nCases]  = generateAllOpts(OptsStruct);

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
            OptsStruct.tol   = {1e-3};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStruct);

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
            OptsStruct.refine      = {0};
            OptsStruct.refTol      = {1e-6};
            OptsStruct.refMaxiter  = {1e2};
            OptsStruct.lse         = {'sparse'};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStruct);

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
            OptsStruct.real        = {true, false};
            OptsStruct.orth        = {false, 'dgks', 'mgs'};
            OptsStruct.reorth      = {'mgs', false, 'qr'};
            OptsStruct.lse         = {'sparse','full','hess','iterative'};
            OptsStruct.dgksTol     = {1e-12};
            OptsStruct.krylov      = {'standardKrylov','cascadedKrylov'};
            OptsStruct.dgksTol     = {1e3};
            OptsStruct.tollse      = {1e-6};
            OptsStruct.solver      = {'cgs','bicgstab','bicg'};
            OptsStruct.verbose     = {true, false};
            OptsStruct.force       = {true, false};

            [AllOptsCell,nCases]  = generateAllOpts(OptsStruct);

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
            OptsStruct.type        = {'SM','LM',1};
            OptsStruct.orth        = {true, false,'qr'};
            OptsStruct.real        = {true, false,'real'};
            OptsStruct.tol         = {1e-6};
            OptsStruct.dominance   = {0,'analyze','2q','3q'};
            OptsStruct.lse         = {'sparse','full'};
            OptsStruct.subspaceW   = {'eigs','1by1'};

            [AllOptsCell,nCases]  = generateAllOpts(OptsStruct);

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
            OptsStruct.trans       = {'T','H'};

            [AllOptsCell,nCases]  = generateAllOpts(OptsStruct);

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
            OptsStruct.real       = {true,false};

            [AllOptsCell,nCases]  = generateAllOpts(OptsStruct);

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
            OptsStruct.type        = {'tbr', 'adi', 'matchDcGain'};
            OptsStruct.redErr      = {'0', 1};
            OptsStruct.hsvTol      = {1e-15};
            OptsStruct.warnOrError = {'warn','error', 0};
            OptsStruct.lse         = {'sparse', 'full','iterative', 'gauss'}; 
            OptsStruct.rctol       = {1e-15};
            OptsStruct.forceOrder  = {true, false};


            [AllOptsCell,nCases]  = generateAllOpts(OptsStruct);

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
            OptsStruct.qm0     = {2*length(s0)};
            OptsStruct.s0m     = {[s0,s0]};
            OptsStruct.maxiter = {15};
            OptsStruct.tol     = {1e-3};
            OptsStruct.stopCrit= {'s0', 'sysr', 'sysm', 'combAny', 'combAll'};

            OptsStruct.verbose         = {true, false};
            OptsStruct.plot            = {true, false};
            OptsStruct.suppressWarn    = {true, false};
            OptsStruct.updateModel     = {'new', 'all'};
            OptsStruct.modelTol        = {1e-2};
            OptsStruct.clearInit       = {true, false};
            OptsStruct.stableModelFct  = {true, false};

            OptsStruct.irka.suppressverbose = {true, false};
            OptsStruct.irka.stopCrit        = {'combAny', 's0', 'sysr', 'combAll'};
            OptsStruct.lse                  = {'sparse', 'full', 'hess'};
            OptsStruct.irka.tol             = {1e-6};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStruct);

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
            OptsStruct.warn            = {true, false};
            OptsStruct.w               = {[]}; % [{''} / '{wmin,wmax}' / vector of frequencies]
            OptsStruct.zeroThres       = {1e-4};
            OptsStruct.cure.verbose    = {true, false};
            OptsStruct.cure.redfun     = {'spark', 'irka', 'rk+pork'};
            OptsStruct.cure.fact       = {'V', 'W'};
            OptsStruct.cure.nk         = {2};
            OptsStruct.cure.stop       = {'normROM', 'nmax', 'h2Error'};
            OptsStruct.cure.stopval    = {1e-6};
            OptsStruct.cure.initMode   = {'zero', 'sm', 'lm', 'slm'};
            OptsStruct.cure.initN      = {5*2};
            OptsStruct.cure.fact       = {'V'};
            
            OptsStruct.cure.test       = {true, false};
            OptsStruct.cure.gif        = {true, false};
            
            OptsStruct.cure.maxIter    = {20};
            OptsStruct.cure.checkEVB   = {true};
            OptsStruct.cure.sEVBTol    = {1e-16};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStruct);

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
            OptsStruct.refine      = {0};
            OptsStruct.refTol      = {1e-6};
            OptsStruct.refMaxiter  = {1e2};
            OptsStruct.lse         = {'sparse'};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStruct);

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
            OptsStruct.maxiter     = {50}; 
            OptsStruct.tol         = {1e-3}; 
            OptsStruct.verbose     = {true, false};
            OptsStruct.stopCrit    = {'combAny', 's0', 'sysr', 'combAll'};
            OptsStruct.suppressverbose = {true, false};
            OptsStruct.lyapchol    = {'', 'adi', 'builtIn'};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStruct);

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
            OptsStruct.updateModel = {'new', 'all', 'new', 'lean'};
            OptsStruct.modelTol    = {1e-2};
            OptsStruct.plot        = {true, false};
            OptsStruct.tol         = {1e-2};

            [AllOptsCell,nCases]  = generateAllOpts(OptsStruct);

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
            OptsStruct.qm0           = {2*length(s0)};
            OptsStruct.s0m           = {[s0,s0]};
            OptsStruct.maxiter       = {15};
            OptsStruct.tol           = {1e-3};
            OptsStruct.stopcrit      = {'combAny', 's0', 'sysr', 'sysm', 'combAny', 'combAll'};
            OptsStruct.verbose       = {true, false};
            OptsStruct.plot          = {true, false};
            OptsStruct.suppressWarn  = {true, false};
            OptsStruct.updateModel   = {'new', 'all'};
            OptsStruct.modelTol      = {1e-2};
            OptsStruct.clearInit     = {true, false};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStruct);

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
            OptsStruct.rk      = {'twoSided', 'input', 'output'}
            OptsStruct.tol     = {1e-2};
            OptsStruct.maxIter = {100};
            OptsStruct.lse     = {'sparse', 'full', 'hess', 'gauss', 'iterative'};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStruct);

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
            OptsStruct.rk = {'twoSided', 'input', 'output'};
            OptsStruct.lse = {'sparse', 'full', 'hess', 'gauss', 'iterative'};
            
            [AllOptsCell,nCases]  = generateAllOpts(OptsStruct);

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
            OptsStruct.zeroThres  = {1e-6}; 

            OptsStruct.spark.type        = {'model', 'standard'};
            OptsStruct.spark.test        = {true, false};
            OptsStruct.spark.verbose     = {true, false};
            OptsStruct.spark.mfe         = {5e3};
            OptsStruct.spark.mi          = {5e3};
            OptsStruct.spark.xTol        = {1e-10};
            OptsStruct.spark.fTol        = {1e-10};
            OptsStruct.spark.modelTol    = {1e-5};
            OptsStruct.spark.pork        = {'V', 'W'};

            OptsStruct.mespark.ritz      = {true, false};
            OptsStruct.mespark.pertIter  = {5};
            OptsStruct.mespark.maxIter   = {20};

            [AllOptsCell,nCases]        = generateAllOpts(OptsStruct);

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