function MORLab_gettingStarted

% MORLAB_GETTINGSTARTED - Introductory demo to MORLab toolbox 
% ------------------------------------------------------------------
% USAGE:  This function can be executed as it is and will guide the 
% user in the command window through an introductory journey in the 
% capabilities and advantages of using sparse state space (sss) 
% objects and the related model reduction functions.
%
% See also SSS, SSS\BODE.
%
% ------------------------------------------------------------------
% This file is part of MORLab, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute 
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                      -> MORLab@tum.de <-
% ------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Last Change:  14 Jul 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%%  Initialization
clear, clc
fprintf('Starting demo execution: MORLab_gettingStarted...\n\n'); 

rule
fprintf(['\t MORLab - A sparse state space and model reduction toolbox\n']); 
rule
fprintf(['\tDeveloped at Chair of Automatic Control, TUM\n']);
fprintf(['\tNote: for academic use only.\n']);
rule

%%  PART I: Sparse State Space (SSS)
fprintf(['\t PART I: Sparse State Space (SSS) objects\n\n']);

fprintf('Large scale linear systems are usually modelled by sparse\n');
fprintf('matrices, i.e. matrices with only a small number of nonzero\n');
fprintf('elements compared to the dimensions. The advantage of using\n');
fprintf('sparse matrices in terms of storage and computational\n');
fprintf('complexity is huge.\n');

try kDouble = 0; while 1, kDouble = kDouble+1; bla = eye(10^kDouble); end     
catch err; end
try kSparse = 0; while 1, kSparse = kSparse+1; bla = speye(10^kSparse); end     
catch err; end
customPause

fprintf('For example, on your computer the maximum size an identity\n');
fprintf('matrix can have if stored as "ful" is 10^%i, wherease\n',kDouble);
fprintf('a sparse identity matrix can reach the dimension 10^%i.\n',kSparse);
customPause

fprintf('Unfortunately, the DynamicSystem class in MATLAB, in particular\n');
fprintf('state space (ss or dss) objects, do not support sparse matrices.\n');
fprintf('As a result, in the large scale setting, where the storing of \n');
fprintf('system matrices as "full" is not possible and/or desirable,\n');
fprintf('the basic system analysis functions such as bode, bodemag, eig \n');
fprintf('freqresp etc., but even model reduction functions as balred, \n');
fprintf('balancmr, modred etc. cannot be used \n')
customPause

fprintf('MORLab includes the definition of sparse state space (sss) objects\n');
fprintf('and the sparsity and large-scale-optimized implementation of some \n');
fprintf('of the most common function for dynamic systems objects.\n\n');
customPause

fprintf('Let us begin with an example that illustrates the capabilities \n');
fprintf('of the MORLab toolbox and in particular sss-objects.\n');

%   *Selection of a benchmark model
result = input(prompt)

function rule
fprintf('--------------------------------------------------------------\n');

function customPause
fprintf('...\n');
pause