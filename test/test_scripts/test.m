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
% Last Change:  04 Sep 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

clear all; 
clc;
import matlab.unittest.TestSuite;

%% Test all unittest-files in folder 'testScripts'
suite = TestSuite.fromFolder('testScripts');


%% Test specific unittest-files
% Available tests:
suite1=TestSuite.fromFile('testArnoldi.m');
suite2=TestSuite.fromFile('testRk.m');
suite3=TestSuite.fromFile('testIrka.m');
suite4=TestSuite.fromFile('testTbr.m'); 
% suite5=TestSuite.fromFile('testScripts\testModal.m'); %not yet

% Add/remove suiteX (e.g. [suite1, suite3] to run testArnoldi and testIRKA)
% suite=[ suite4];


%% Run and show results
result = run(suite);
disp(result);
