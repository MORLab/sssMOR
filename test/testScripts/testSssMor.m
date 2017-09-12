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
%                   -> morlab@rt.mw.tum.de <-
%------------------------------------------------------------------
% Authors:      Lisa Jeschek, Jorge Luiz Moreira Silva
% Last Change:  11 Feb 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

import matlab.unittest.TestSuite;

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
testPathSssMor = p(1:k(end)-1);

% load benchmarks
if Opts.loadBench == 1
    testPathSss=loadBenchmarks(Opts);
end

%% Test specific unittest-files
% Classic MOR:
suiteClassic = [TestSuite.fromFile(fullfile(testPathSssMor,'testArnoldi.m')),...
                TestSuite.fromFile(fullfile(testPathSssMor,'testDemoSssMor.m')),...
                TestSuite.fromFile(fullfile(testPathSssMor,'testIrka.m')),...
                TestSuite.fromFile(fullfile(testPathSssMor,'testIsH2opt.m')),...
                TestSuite.fromFile(fullfile(testPathSssMor,'testIsrk.m')),...
                TestSuite.fromFile(fullfile(testPathSssMor,'testModal.m')),... 
                TestSuite.fromFile(fullfile(testPathSssMor,'testMoments.m')),...
                TestSuite.fromFile(fullfile(testPathSssMor,'testProjectiveMor.m')),...
                TestSuite.fromFile(fullfile(testPathSssMor,'testRk.m')),...
                TestSuite.fromFile(fullfile(testPathSssMor,'testRkOp.m')),...
                TestSuite.fromFile(fullfile(testPathSssMor,'testRkIcop.m')),...
                TestSuite.fromFile(fullfile(testPathSssMor,'testTbr.m'))]; 


% State of the art MOR:
suiteStateOfTheArt = [TestSuite.fromFile(fullfile(testPathSssMor,'testSylvester.m')),...
                      TestSuite.fromFile(fullfile(testPathSssMor,'testPork.m')),...
                      TestSuite.fromFile(fullfile(testPathSssMor,'testSpark.m')),...
                      TestSuite.fromFile(fullfile(testPathSssMor,'testCure.m')),...
                      TestSuite.fromFile(fullfile(testPathSssMor,'testCirka.m')),...
                      TestSuite.fromFile(fullfile(testPathSssMor,'testModelFctMor.m'))];
                  
% ssRed class
suiteSsRed = [TestSuite.fromFile(fullfile(testPathSssMor,'testSsRed.m')),...
              TestSuite.fromFile(fullfile(testPathSssMor,'testStabsep.m')),...
              TestSuite.fromFile(fullfile(testPathSssMor,'testL2norm.m'))  ];


suiteAll=[suiteClassic, suiteStateOfTheArt, suiteSsRed];

%% Run and show results
% Choose between suiteClassic, suiteStateOfTheArt, suiteAll
eval(sprintf('result = run(%s);',Opts.suite));
disp(result);

if Opts.loadBench==1
	(fullfile(testPathSss,'benchmarksSysCell.mat'));
end
end

