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
        function runInitializeShigts(testCase)

        end
        function runIsH2opt(testCase)

        end
        function runMoments(testCase)

        end
        function runSsRed(testCase)

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
                    %  test system
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
                    %  test system
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
                    %  test system

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

            h = waitbar(0,'projectiveMor: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOptsCell{kOpts};
                for i=1:length(testCase.sysCell)
                    %  test system
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
                    %  test system
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

        end
        function runCure(testCase)

        end
        function runIrka(testCase)

        end
        function runIsrk(testCase)

        end
        function runModelFct(testCase)

        end
        function runModelFctMor(testCase)

        end
        function runRkIcop(testCase)

        end
        function runRkOp(testCase)

        end
        function runSpark(testCase)

        end
    end
end