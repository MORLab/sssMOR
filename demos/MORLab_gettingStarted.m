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
fprintf('of the most common function for dynamic systems objects.\n');
customPause


fprintf('\nLet us begin with an example that illustrates the capabilities \n');
fprintf('of the MORLab toolbox and in particular sss-objects.\n');

%   *Selection of a benchmark model
fprintf(['Choose one of the following benchmark models by pressing...\n',...
                '\t "1" \t ..for the Los Angeles Hospital model (N=48)\n',...
                '\t "2" \t ..for a 1D linear beam model (N=348)\n']);
modelChosen = input('Model choice: ');

if isempty(modelChosen),modelChosen = 0;end
isValid = 0;
while ~isValid
    switch modelChosen
        case 1 %build
            sysName = 'build';
            load ../benchmarks/build.mat
            E = speye(size(A));
            D = 0;
            isValid = 1;
        case 2 %beam
            sysName = 'beam';
            load ../benchmarks/beam.mat
            E = speye(size(A));
            D = 0;
            isValid = 0;
        otherwise    
            warning('Your input was not valid')
            fprintf(['\nChoose one of the following benchmark models by pressing...\n',...
                '\t "1" \t ..for the Los Angeles Hospital model (N=48)\n',...
                '\t "2" \t ..for a 1D-beam model (N=348)\n']);
            modelChosen = input('Model choice: ');
    end
end

%   *Display some information on the model
fprintf('\nYou chose the "%s" model. \n',sysName);
fprintf('The order of the model is %i. \n',size(A,1));
fprintf('The model has %i input(s) and %i output(s). \n',size(B,2), size(C,1));
customPause


%   *Create an sss-object sys
fprintf('\nThe dynamics of the system are described by the implicit\n');
fprintf('state space realization with sparse matrices:\n');
fprintf('\tE * d/dt(x) = A x + B u\n');
fprintf('\t         y  = C x + D u\n');
customPause

fprintf('Using MORLab, the dynamic system can be stored as sss-object calling:\n');
fprintf('>> sys = sss(A,B,C,D,E)\n');
sys = sss(A,B,C,D,E);
sysDss = dss(full(A),full(B),full(C),full(D),full(E));

W = whos('sys','sysDss');
storagePerc = W(1).bytes/W(2).bytes*100;
fprintf('(Note that the storage required for the sss-object is %03.2f%% \n',storagePerc);
fprintf('of the one needed for the respective ss-object.)\n');
customPause

fprintf('The sparsity pattern of E and A can be plotted by calling\n');
fprintf('>> spysys(sys)\n');
fprintf('\t\tor\n');
fprintf('>> spysys(sys.E, sys.A)\n');
spysys(sys);
customPause

%*  Bode plot
fprintf('Another function that can be useful to analyze the dynamic\n');
fprintf('of the given system is the "bode" function, called using:\n');
fprintf('>> bode(sys)\n');
figure; tic;bode(sys);tBodeSss = toc;

% fh = figure; tic, bode(sysDss); tBodeSs = toc;close(fh)
% tBodePerc = tBodeSss/tBodeSs*100;
% fprintf('(Note that the time required to plot the sss-object is %03.2f%% \n',tBodePerc);
% fprintf('of the one needed for the respective ss-object.)\n');
customPause

%*  System norms
%H2
fprintf('System norms can also be computed using the "norm" function.\n');
fprintf('H2-norm:\n');
fprintf('>> norm(sys)\n');
fprintf('\tor\n');
fprintf('>> norm(sys,2)\n');

    tic;H2norm = norm(sys);tNorm2Sss = toc;
fprintf('The H2-norm of the chosen system is: %e\n',H2norm);

    tic, norm(sysDss); tNorm2Ss = toc;
tH2Perc = tNorm2Sss/tNorm2Ss*100;
fprintf('(Note that the time required to compute the norm of the');
fprintf('sss-object is %03.2f%% of the one needed for the\n',tH2Perc);
fprintf('respective ss-object.)\n');

%Hinf
fprintf('System norms can also be computed using the "norm" function.\n');
fprintf('H2-norm:\n');
fprintf('>> norm(sys)\n');
fprintf('\tor\n');
fprintf('>> norm(sys,2)\n');

    tic;H2norm = norm(sys);tNorm2Sss = toc;
fprintf('The H2-norm of the chosen system is: %e\n',H2norm);

    tic, norm(sysDss); tNorm2Ss = toc;
tH2Perc = tNorm2Sss/tNorm2Ss*100;
fprintf('(Note that the time required to compute the norm of the');
fprintf('sss-object is %03.2f%% of the one needed for the\n',tH2Perc);
fprintf('respective ss-object.)\n');


customPause

%*  Time domain analysis.
fprintf('System norms can also be computed using the "norm" function.\n');
fprintf('H2-norm:\n');
fprintf('>> norm(sys)\n');
fprintf('\tor\n');
fprintf('>> norm(sys,2)\n');

    tic;H2norm = norm(sys);tNorm2Sss = toc;
fprintf('The H2-norm of the chosen system is: %e\n',H2norm);

    tic, norm(sysDss); tNorm2Ss = toc;
tH2Perc = tNorm2Sss/tNorm2Ss*100;
fprintf('(Note that the time required to compute the norm of the');
fprintf('sss-object is %03.2f%% of the one needed for the\n',tH2Perc);
fprintf('respective ss-object.)\n');
customPause








function rule
fprintf('--------------------------------------------------------------\n');

function customPause
fprintf('...\n');
pause