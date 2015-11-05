% TEST - testing of sssMOR functionality
%
% Description:
%   This file provides test functions for the sssMOR toolbox. All tests
%   contained in the folder 'testScripts' can be executed or the tests can
%   be run seperately by choosing the test-suites.
%
%   Required systems: build, beam, fom ,eady, random, LF10 (IMTEK), 
%                     SpiralInductorPeec (IMTEK)
% ------------------------------------------------------------------
%   This file is part of sssMOR, a Sparse State Space, Model Order
%   Reduction and System Analysis Toolbox developed at the Institute 
%   of Automatic Control, Technische Universitaet Muenchen.
%   For updates and further information please visit www.rt.mw.tum.de
%   For any suggestions, submission and/or bug reports, mail us at
%                     -> sssMOR@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
%               Lisa Jeschek
% Last Change:  07 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

clear all; 
clc;
import matlab.unittest.TestSuite;

%%  Change to testScripts folder
testCase.Path = pwd; %original
            
p = mfilename('fullpath'); k = strfind(p, '\test'); 
testpath = p(1:k(end)-1);
cd(testpath);
%%
% 
% * ITEM1
% * ITEM2
% 

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
suite8=TestSuite.fromFile('testMomentsAll.m'); %testMoments using all benchmarks in a local folder
suite9=TestSuite.fromFile('testMtimes.m');
suite10=TestSuite.fromFile('testNorm.m');
suite11=TestSuite.fromFile('testPzmap.m');
suite12=TestSuite.fromFile('testPlus.m');

% suite11=TestSuite.fromFile('testSSS.m'); %not yet

% Add/remove suiteX (e.g. [suite1, suite3] to run testArnoldi and testIRKA)
suite=[suite1,suite2,suite3,suite4,suite5,suite6,suite7,suite8,suite9,suite10,suite11,suite12];


%% Run and show results
result = run(suite);
disp(result);

%% Go back to original folder
cd(testCase.Path);
