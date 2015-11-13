function result = test(Opts)
% TEST - testing of sssMOR functionality
%
% Description:
%   This file provides test functions for the sssMOR toolbox. All tests
%   contained in the folder 'testScripts' can be executed or the tests can
%   be run seperately by choosing the test-suites.
%
%   Required systems: build, beam, fom, eady, random, LF10 (IMTEK), 
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
% Last Change:  10 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

clear;
clc;
import matlab.unittest.TestSuite;

%%  Change to testScripts folder
testCase.Path = pwd; %original
            
p = mfilename('fullpath'); k = strfind(p, '\test'); 
testpath = p(1:k(end)-1);
cd(testpath);

%% Choose benchmarks
% Default benchmarks
Def.option = 'light'; %'light','full','heavy'
Def.size = 400'; %'light': only benchmarks with sys.n<=Opts.size are tested
                %'heavy': only benchmarks with sys.n>Opts.size are tested
Def.number = 3; %choose maximum number of tested benchmarks

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

%% Load benchmarks
%the directory "benchmark" is in sssMOR
p = mfilename('fullpath'); k = strfind(p, 'test\'); 
pathBenchmarks = [p(1:k-1),'benchmarks'];
cd(pathBenchmarks);
badBenchmarks = {'LF10.mat','beam.mat','random.mat',...
    'SpiralInductorPeec.mat'};  

% load files
files = dir('*.mat'); 
benchmarksSysCell=cell(1,Opts.number);
nLoaded=1; %count of loaded benchmarks
disp('Loaded systems:');

warning('off');
for i=1:length(files)
    if nLoaded<Opts.number+1
        switch(Opts.option)
            case 'light'
                sys = loadSss(files(i).name);
                if ~any(strcmp(files(i).name,badBenchmarks)) && size(sys.A,1)<=Opts.size
                    benchmarksSysCell{nLoaded}=sys;
                    nLoaded=nLoaded+1;
                    disp(files(i).name);
                end
            case 'full' 
                benchmarksSysCell{nLoaded} = loadSss(files(i).name);
                nLoaded=nLoaded+1;
                disp(files(i).name);
            case 'heavy'
                sys = loadSss(files(i).name);
                if any(strcmp(files(i).name,badBenchmarks)) || size(sys.A,1)>Opts.size
                    benchmarksSysCell{nLoaded}=sys;
                    nLoaded=nLoaded+1;
                    disp(files(i).name);
                end
            otherwise
                error('Benchmark option is wrong.');
        end
    end
end
benchmarksSysCell(nLoaded:end)=[];
warning('on');

% change path back and save loaded systems
cd(testpath);
save('benchmarksSysCell.mat');

%% Test all unittest-files in current folder
% suite = TestSuite.fromFolder(pwd);


%% Test specific unittest-files
% Available tests:
suite1=TestSuite.fromFile('testArnoldi.m');
suite2=TestSuite.fromFile('testRk.m');
suite3=TestSuite.fromFile('testIrka.m');
suite4=TestSuite.fromFile('testTbr.m'); 
suite5=TestSuite.fromFile('testModal.m'); 
suite6=TestSuite.fromFile('testMoments.m'); 
suite7=TestSuite.fromFile('testIsH2opt.m');
suite8=TestSuite.fromFile('testMtimes.m');
suite9=TestSuite.fromFile('testNorm.m');
suite10=TestSuite.fromFile('testPzmap.m');
suite11=TestSuite.fromFile('testPlus.m');
suite12=TestSuite.fromFile('testResidue.m');
suite13=TestSuite.fromFile('testSigma.m');
suite14=TestSuite.fromFile('testStep.m');
suite15=TestSuite.fromFile('testSs.m');
suite16=TestSuite.fromFile('testAppend.m');
suite17=TestSuite.fromFile('testBode.m');
suite18=TestSuite.fromFile('testEig.m');
suite19=TestSuite.fromFile('testEigs.m');
suite20=TestSuite.fromFile('testFreqresp.m');
suite21=TestSuite.fromFile('testImpulse.m');
suite22=TestSuite.fromFile('testIsstable.m');
suite23=TestSuite.fromFile('testDiag.m');
suite24=TestSuite.fromFile('testDecayTime.m');
suite25=TestSuite.fromFile('testIssd.m');
% suite26=TestSuite.fromFile('testSim.m');
suite27=TestSuite.fromFile('testSylvester.m');
suite28=TestSuite.fromFile('testPork.m');
suite29=TestSuite.fromFile('testSpark.m');
suite30=TestSuite.fromFile('testCure.m');
suite31=TestSuite.fromFile('testConnectSss.m');


% Add/remove suiteX (e.g. [suite1, suite3] to run testArnoldi and testIRKA)

% suite=[suite1,suite2,suite3,suite4,suite5,suite6,suite7,suite8,suite9,suite10,...
% suite11,suite12,suite13,suite14,suite15,suite16, suite17, suite18, suite19,...
% suite20, suite21, suite22, suite23];

suite=[suite1,suite2,suite3,suite4,suite5,suite6,suite7,suite8,suite9,suite10,...
suite11,suite12,suite13,suite14,suite15,suite16, suite17,...
suite20, suite21, suite22, suite25, suite27, suite28, suite29, suite30];

%% Run and show results
result = run(suite);
disp(result);
delete('benchmarksSysCell.mat');

%% Go back to original folder
cd(testCase.Path);
