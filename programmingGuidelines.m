% PROGRAMMING GUIDELINES
%
% USAGE:  This script summarizes the naming and coding conventions 
% used for the development of the MORLab Toolbox. All developers must
% abide by these guidelines.
%
% REFERENCES:
% [1] Johnson (2014): MATLAB Style Guidelines 2.0
% [2] Hung (1999): The Pragmatic Programmer
% 
% ------------------------------------------------------------------
% This file is part of sssMOR, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute 
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                    -> morlab@rt.mw.tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Last Change:  26 April 2015
% Copyright 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%%  1) Naming conventions
%   general variable and function names should be mixed-case, starting with
%   a lower case, and be self explainatory. Acronyms should be lower case.
%   Underscores should be avoided! Abbreviations should be avoided unless
%   they are commonly used in the domain
krylovSubspace, shift, reducedOrderModel, html, tum
rom = reducedOrderModel %examples of acceptable abbreviations
fom = fullOrderModel

%   amounts, number of objects: n<Objects>
nShifts, nPoints

%   avoid plurals
shift % one
shiftArray, shiftVec % many

%   indices: i<Object>
iShift, iFile, jPosition
for iPoint = 1:nPoints %example of nested loops
    for jFile = 1:nFiles
    end
end 

%   booleans: use positive bolean names
isFound %instead of isNotFound

%   constants (that should not be changed) are upper case. They can be 
%   preceded by a prefix that specifies their group
MAX_ITERATIONS, COLOR_RED, COLOR_BLUE

%   Structures begin with upper case letter. Fieldnames do not repeat the
%   name of the structure
Segment.length, Options.MAX_ITERATIONS, GeneralOptions.fieldName

%   Functions that return only one value are named by that value
mean, moment, maxError

%%  2) Definition of execution options
%   Whenever a function can be executed with different parameters or even
%   in different modes, the execution options should be passed at the end
%   of the inputs such as
      sysr = cure(sys,..,Opts)
%   Opts is a structure containing some of the parameters that the function
%   accepts.
%   Default paramters are defined at the beginning of the function in which
%   they are used with a default structure
      Def.<fieldname> = ...
%   Then, the input should be parsed and the current options structure
%   updated with the function parseOpts:
      if ~exist('Opts','var') || isempty(Opts)
            Opts = Def;
      else
            Opts = parseOpts(Opts,Def);
      end       
      
%%  3) Definition of empty (sparse) matrices for memory allocation
A = sparse(n,n); %zeros(n) can go out of memory

