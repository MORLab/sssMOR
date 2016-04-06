function result = testSssMor(Opts)
% TEST - testing of sssMOR functionality
%
% Description:
%   This file provides test functions for the sssMOR toolbox. All tests
%   contained in the folder 'testScripts' can be executed or the tests can
%   be run seperately by choosing the test-suites.
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

import matlab.unittest.TestSuite;

%%  Change to testScripts folder
testCase.Path = pwd; %original

%% Choose benchmarks
% Default benchmarks
Def.cond = 'good'; % condition of benchmarks: 'good','bad','all'
                   % 'bad': LF10, beam, random, SpiralInductorPeec
                   % 'good': all benchmarks that are not 'bad'
Def.minSize     = 0; % test benchmarks with sys.n >= minSize
Def.maxSize     = 400; % test benchmarks with sys.n <= minSize
Def.number      = 3; % choose maximum number of tested benchmarks
Def.loadBench   = 1; %loadBenchmarks (for testAll)
Def.suite       = 'suiteAll';

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

% sssMor testScripts
p = mfilename('fullpath'); k = strfind(p, fullfile(filesep,'test')); 
testpathSssMor = p(1:k(end)-1);

% load benchmarks and change to folder 'testScripts'
if Opts.loadBench == 1
    testpathSss=loadBenchmarks(Opts);
end

cd(testpathSssMor);

%% Test specific unittest-files
% Classic MOR:
suiteClassic = [TestSuite.fromFile('testArnoldi.m'),...
                TestSuite.fromFile('testRk.m'),...
                TestSuite.fromFile('testIrka.m'),...
                TestSuite.fromFile('testTbr.m'), ...
                TestSuite.fromFile('testModal.m'),... 
                TestSuite.fromFile('testMoments.m'),... 
                TestSuite.fromFile('testMoments.m')]; 


% State of the art MOR:
suiteStateOfTheArt = [TestSuite.fromFile('testSylvester.m'),...
                      TestSuite.fromFile('testPork.m'),...
                      TestSuite.fromFile('testSpark.m'),...
                      TestSuite.fromFile('testCure.m')];

suiteAll=[suiteClassic, suiteStateOfTheArt];

%% Run and show results
% Choose between suiteClassic, suiteStateOfTheArt, suiteAll
eval(sprintf('result = run(%s);',Opts.suite));
disp(result);

if Opts.loadBench==1
    cd(testpathSss);
	('benchmarksSysCell.mat');
end

%% Go back to original folder
cd(testCase.Path);
end

