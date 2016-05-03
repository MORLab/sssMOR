function result = test(Opts)
% TEST - testing of sssMOR and sss functionality
%
% Description:
%   This file provides test functions for the sssMOR toolbox. All tests
%   contained in the folders 'testScripts' of sss and sssMor can be 
%   executed.
%
%   Required systems: building, beam, fom, eady, random, LF10 (IMTEK), 
%                     SpiralInductorPeec (IMTEK), rail_1357 (IMTEK)
%------------------------------------------------------------------
% This file is part of sssMOR, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute 
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                   -> sssMOR@rt.mw.tum.de <-
%------------------------------------------------------------------
% Authors:      Lisa Jeschek, Jorge Luiz Moreira Silva
% Last Change:  11 Feb 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

% Default benchmarks
Def.cond = 'good'; % condition of benchmarks: 'good','bad','all'
                   % 'bad': LF10, beam, random, SpiralInductorPeec
                   % 'good': all benchmarks that are not 'bad'
Def.minSize = 0; % test benchmarks with sys.n >= minSize
Def.maxSize = 400; % test benchmarks with sys.n <= minSize
Def.number = 3; % choose maximum number of tested benchmarks
Def.loadBench = 0;

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
    Opts.loadBench = 0;
end

testPath=loadBenchmarks(Opts);

% Run testSssMor
disp([10,'Testing SSSMOR...']);
resultSssMor=testSssMor(Opts);
disp([10,'Testing SSS...']);
resultSss=testSss(Opts);

disp([10,'Result of testing SSSMOR:']);
disp(resultSssMor)
disp(table(resultSssMor));
disp([10,'Result of testing SSS:']);
disp(resultSss);
disp(table(resultSss));
result={resultSssMor, resultSss};

delete(fullfile(testPath,'benchmarksSysCell.mat'));
end

