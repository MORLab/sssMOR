function  sss_gettingStarted
% SSS_GETTINGSTARTED - Introductory demo to sssMOR toolbox 
% 
% Syntax:
%       SSS_GETTINGSTARTED
% 
% Description:
%       This function can be executed as it is and will guide the 
%       user in the command window through an introductory journey in the 
%       capabilities and advantages of using sparse state space (sss) 
%       objects.
%
% See Also:
%       MOR_gettingStarted
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sss">sss</a>, a Sparse State-Space and System Analysis 
% Toolbox developed at the Chair of Automatic Control in collaboration
% with the Chair of Thermofluid Dynamics, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  23 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Initialization
clear, clc
warning('off','all');
fprintf('Starting demo execution: sssMOR_gettingStarted...\n\n'); 

rule
fprintf(['\t sssMOR - A sparse state space and model reduction toolbox\n']); 
rule
fprintf(['\tDeveloped at Chair of Automatic Control, TUM\n']);
fprintf(['\tNote: for academic use only.\n']);
rule

%%  Sparse State Space (SSS)
fprintf(['\t Sparse State Space (SSS) objects\n\n']);

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
fprintf('matrix can have if stored as...\n');
fprintf('\t... "full" is 10^%i\n',kDouble);
fprintf('\t... sparse is 10^%i.\n',kSparse);
customPause

fprintf('Unfortunately, the DynamicSystem class in MATLAB, in particular\n');
fprintf('state space (ss or dss) objects, do not support sparse matrices.\n');
fprintf('As a result, in the large scale setting, where the storing of \n');
fprintf('system matrices as "full" is not possible and/or desirable,\n');
fprintf('the basic system analysis functions such as bode, bodemag, eig \n');
fprintf('freqresp etc., but even model reduction functions as balred, \n');
fprintf('balancmr, modred etc. cannot be used \n')
customPause

fprintf('sssMOR includes the definition of sparse state space (sss) objects\n');
fprintf('and the sparsity and large-scale-optimized implementation of some \n');
fprintf('of the most common function for dynamic systems objects.\n');
customPause


fprintf('\nLet us begin with an example that illustrates the capabilities \n');
fprintf('of the sssMOR toolbox and in particular sss-objects.\n');

%   *Selection of a benchmark model
[sysName,A,B,C,D,E] = selectModel;

%   *Display some information on the model
fprintf('\nYou chose the "%s" model. \n',sysName);
% fprintf('The order of the model is %i. \n',size(A,1));
% fprintf('The model has %i input(s) and %i output(s). \n',size(B,2), size(C,1));
% customPause


%   *Create an sss-object sys
fprintf('\nThe dynamics of the system are described by the explicit\n');
fprintf('state space realization with sparse matrices:\n');
fprintf('\t    d/dt(x) = A x + B u\n');
fprintf('\t         y  = C x + D u\n');
customPause

fprintf('Using sssMOR, the dynamic system can be stored as sss-object calling:\n');
fprintf('>> sys = sss(A,B,C,D,E)\n');
sys = sss(A,B,C,D,E);
sysDss = dss(full(A),full(B),full(C),full(D),full(E));
clear A B C D E

W = whos('sys','sysDss');
storagePerc = W(1).bytes/W(2).bytes*100;
fprintf('(Note that the storage required for the sss-object is %03.2f%% \n',storagePerc);
fprintf('of the one needed for the respective ss-object.)\n');
customPause

fprintf('To display some information about the model, use:\n');
fprintf('>> disp(sys)\n');
disp(sys)
customPause

fprintf('The sparsity pattern of E and A can be plotted by calling\n');
fprintf('>> spy(sys)\n');
fprintf('\t\tor\n');
fprintf('>> spy(sys.E, sys.A)\n');
spy(sys);
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
% TODO: Why is our  bode so mach slower than the builtin?
customPause

%*  System norms
%H2
fprintf('System norms can also be computed using the "norm" function.\n');
fprintf('H2-norm:\n');
fprintf('\t>> norm(sys)\n');
fprintf('\t\tor\n');
fprintf('\t>> norm(sys,2)\n');

    tic;H2norm = norm(sys);tNorm2Sss = toc;
fprintf('\tThe H2-norm of the chosen system is: %e\n',H2norm);

    tic, norm(sysDss); tNorm2Ss = toc;
tH2Perc = tNorm2Sss/tNorm2Ss*100;
fprintf('\t(Note that the time required to compute the norm of the \n');
fprintf('\tsss-object is %03.2f%% of the one needed for the\n',tH2Perc);
fprintf('\trespective ss-object.)\n');

% %Hinf
% fprintf('HInf-norm:\n');
% fprintf('\t>> norm(sys,Inf)\n');
% 
%     tic;HInfnorm = norm(sys,Inf);tNormInfSss = toc;
% fprintf('\tThe HInf-norm of the chosen system is: %e\n',HInfnorm);
% 
%     tic, norm(sysDss,'Inf'); tNormInfSs = toc;
% tHInfPerc = tNormInfSss/tNormInfSs*100;
% fprintf('\t(Note that the time required to compute the norm of the\n');
% fprintf('\tsss-object is %03.2f%% of the one needed for the\n',tHInfPerc);
% fprintf('\trespective ss-object.)\n');
% customPause

%*  Time domain analysis.
fprintf('Often times we wish to analyze the dynamic response of the \n');
fprintf('system in time domain. For this purpose, the "step" function\n');
fprintf('can be used by calling:\n');
fprintf('>> step(sys)\n');
figure;tic; step(sys);tStepSss = toc;

    fh = figure; tic, step(sysDss); tStepSs = toc;close(fh)
tStepPerc = tStepSss/tStepSs*100;
fprintf('(Note that the time required to plot the sss-object is %03.2f%% \n',tStepPerc);
fprintf('of the one needed for the respective ss-object.)\n');
customPause

rule
fprintf('This was just a short introduction to sss-objects and their usage.\n');
fprintf('For more details on sss-objects and the functions implemented,\n');
fprintf('please refer to the documentation of the sss-class.\n');
rule
demoEnd
close all;
rule
function rule
fprintf('--------------------------------------------------------------\n');
function customPause
fprintf('\n');
pause
function [sysName,A,B,C,D,E] = selectModel
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
            LoadData = load('build.mat'); %assumes it is in path
            
            A = LoadData.A; B = LoadData.B; C = LoadData.C;
            E = speye(size(A));
            D = 0;
            isValid = 1;
        case 2 %beam
            sysName = 'beam';
            LoadData = load('beam.mat');%assumes it is in path
            
            A = LoadData.A; B = LoadData.B; C = LoadData.C;
            E = speye(size(A));
            D = 0;
            isValid = 1;
        otherwise    
            warning('on','all');
            warning('Your input was not valid')
            warning('off','all');
            fprintf(['\nChoose one of the following benchmark models by pressing...\n',...
                '\t "1" \t ..for the Los Angeles Hospital model (N=48)\n',...
                '\t "2" \t ..for a 1D-beam model (N=348)\n']);
            modelChosen = input('Model choice: ');
            if isempty(modelChosen),modelChosen = 0;end
    end
end
