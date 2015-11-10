% UNITTEST TEMPLATE
%
% Description:
%   Tests are to be created using the Matlab Unit Test Framework. 
%   Tests are defined as classes. This is an example
%
%   Introductory video available at
%   url = 'http://de.mathworks.com/support/2015a/matlab/8.5/demos/matlab-unit-test-framework-in-release-2013a.html';
%   web(url,'-browser')
% 
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
% Last Change:  04 Oct 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%% Test template

%Definition of a unittest class:
classdef testName < matlab.unittest.TestCase
    
    %Properties (optional)
    properties
          OriginalPath 
    end
    
    %Constructor (optional)
    methods(TestMethodSetup)
        function getBenchmarks(testCase)
            load('benchmarksSysCell.mat');
            testCase.sysCell=benchmarksSysCell;
        end
    end
    
    %Destructor (optional)
    methods(TestMethodTeardown)
    end
    
    %Test functions
    methods (Test)  
        
        function testfunction1(testCase)
            %load benchmark
            load('build.mat');
            s0=5;
            q=2;
            
            %actual solution
            actM = moments(sss(A,B,C), s0, q);
            
            actSolution = actM; %multiple solutions -> cell: e.g {V,W,sysr.A}
            
            %expcected solution
            expM=zeros(1,q);
            temp=(A-s0*eye(size(A))\B);
            expM(1)=C*temp;
            temp=(A-s0*eye(size(A))\temp);
            expM(2)=C*temp;
            
            expSolution = expM;
            
            %compare actual to expected solution
            verifyEqual(testCase, actSolution, expSolution,'RelTol',0.1,'AbsTol',0.0001,...
                'Difference between actual and expected exceeds relative or absolute tolerance');
                %all available test commands can be found at:
                %http://de.mathworks.com/help/matlab/ref/matlab.unittest.qualifications.verifiable-class.html

                %add RelTol (relative tolerance) and/or AbsTol (absolute
                %tolerance) and error message
        end
        
        function testfunction2 (testCase)       
        %add more testfunctions...
        end
    end
end


%% Run test

% Run a single test 
result=run(testName,'testfunction1');
disp(result);

% Run several tests 
import matlab.unittest.TestSuite;
suite1=TestSuite.fromFile('testName1.m');
suite2=TestSuite.fromFile('testName2.m');
suite3=TestSuite.fromFile('testName3.m');
suite=[suite1, suite2, suite3];
result=run(suite);
disp(result);

% Test all unittest-files in current folder(pwd) or path ('testScripts')
import matlab.unittest.TestSuite;
suite = TestSuite.fromFolder(pwd);
result = run(suite);
disp(result);