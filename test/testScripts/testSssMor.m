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
Def.option = 'light'; %'light','full','heavy'
Def.size =400'; %'light': only benchmarks with sys.n<=Opts.size are tested
                %'heavy': only benchmarks with sys.n>Opts.size are tested
Def.number = 3; %choose maximum number of tested benchmarks
Def.loadBench = 1; %loadBenchmarks (for testAll)

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
suite1=TestSuite.fromFile('testArnoldi.m');
suite2=TestSuite.fromFile('testRk.m');
suite3=TestSuite.fromFile('testIrka.m');
suite4=TestSuite.fromFile('testTbr.m'); 
suite5=TestSuite.fromFile('testModal.m'); 
suite6=TestSuite.fromFile('testMoments.m'); 

% State of the art MOR:
suite7=TestSuite.fromFile('testSylvester.m');
suite8=TestSuite.fromFile('testPork.m');
suite9=TestSuite.fromFile('testSpark.m');
suite10=TestSuite.fromFile('testCure.m');

suiteClassic=[suite1,suite2,suite3,suite4,suite5,suite6];
suiteStateOfTheArt=[suite7,suite8,suite9,suite10];
suiteAll=[suite1,suite2,suite3,suite4,suite5,suite6,suite7,suite8,suite9,suite10];

%% Run and show results
% Choose between suiteClassic, suiteStateOfTheArt, suiteAll
result = run(suiteAll);
disp(result);

if Opts.loadBench==1
    cd(testpathSss);
	('benchmarksSysCell.mat');
end

%% Go back to original folder
cd(testCase.Path);
end

