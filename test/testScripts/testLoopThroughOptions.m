classdef testLoopThroughOptions < sssTest
% testModal - testing of modalMor.m
%
% Description:
%  
% Authors:      
% Last Change:  24 Oct 2017
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------
 
    methods(Test)
        function loopThroughOptions(testCase) 
            %run modal mor for all loaded models with all options to check
            %that all combinations are valid
            
            % Define possible opts combinations
            AllOpts.type        = {'SM','LM',1};
            AllOpts.orth        = {true, false,'qr'};
            AllOpts.real        = {true, false,'real'};
            AllOpts.tol         = {1e-6};
            AllOpts.dominance   = {0,'analyze','2q','3q'};
            AllOpts.lse         = {'sparse','full'};
            AllOpts.subspaceW   = {'eigs','1by1'};
            
            [AllOpts,nCases]    = generateAllOpts(AllOpts);
            
            h = waitbar(0,'modalMor: testing all combinations for Opts...');
            try
            for kOpts = 1:nCases
                waitbar(kOpts/nCases,h);
                Opts = AllOpts{kOpts};
                for i=1:length(testCase.sysCell)
                    %  test system
                    sys     = testCase.sysCell{i};
                    if strcmp(Opts.dominance,'2q') && strcmp(testCase.sysCell{i}.Name,'iss')
                        % Option dominance = '2q' produces a numerical error
                        % for the benchmark "iss", because only 2 of the 
                        % requested 8 eigenvalues converge.
                        
                        warning('off')
                        try
                            sysr    = modalMor(sys,4, Opts);
                        end
                        warning('on')
                    else
                        sysr    = modalMor(sys,4, Opts);
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
    end
end