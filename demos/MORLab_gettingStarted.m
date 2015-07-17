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
%                    -> MORLab@rt.mw.tum.de <-
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

fprintf('MORLab includes the definition of sparse state space (sss) objects\n');
fprintf('and the sparsity and large-scale-optimized implementation of some \n');
fprintf('of the most common function for dynamic systems objects.\n');
customPause


fprintf('\nLet us begin with an example that illustrates the capabilities \n');
fprintf('of the MORLab toolbox and in particular sss-objects.\n');

%   *Selection of a benchmark model
[sysName,A,B,C,D,E] = selectModel;

%   *Display some information on the model
fprintf('\nYou chose the "%s" model. \n',sysName);
% fprintf('The order of the model is %i. \n',size(A,1));
% fprintf('The model has %i input(s) and %i output(s). \n',size(B,2), size(C,1));
% customPause


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

%Hinf
fprintf('HInf-norm:\n');
fprintf('\t>> norm(sys,Inf)\n');

    tic;HInfnorm = norm(sys,Inf);tNormInfSss = toc;
fprintf('\tThe HInf-norm of the chosen system is: %e\n',HInfnorm);

    tic, norm(sysDss,'Inf'); tNormInfSs = toc;
tHInfPerc = tNormInfSss/tNormInfSs*100;
fprintf('\t(Note that the time required to compute the norm of the\n');
fprintf('\tsss-object is %03.2f%% of the one needed for the\n',tHInfPerc);
fprintf('\trespective ss-object.)\n');
customPause

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
customPause

%%  PART II: Model Order Reduction
clear

fprintf(['\t PART II: Model Order Reduction (MOR)\n\n']);

fprintf('Besides allowing the definition of sparse dynamic systems, the \n');
fprintf('main goal of the MORLab toolbox is to provide basic and state-of- \n');
fprintf('the-art MOR functions for academic research and teaching.\n');

fprintf('The functions are implemented to work in combination with the sss- \n');
fprintf('class and comprise classical MOR techniques, such as modal reduction, \n');
fprintf('balanced truncation or rational Krylov, as well as more advanced \n');
fprintf('techniques developed in recent years. \n');
customPause

fprintf('\nWe will now demonstrate the main MOR capabilities of the MORLab \n');
fprintf('toolbox through an example.\n');

%   *Selection of a benchmark model
[sysName,A,B,C,D,E] = selectModel;

%   *Display some information on the model
fprintf('\nYou chose the "%s" model. \n',sysName);

fprintf('\nDefine a sss-object and show some information about it.\n');
fprintf('>> sys = sss(A,B,C,D,E)\n');
fprintf('>> disp(sys)\n');
fprintf('>> spysys(sys)\n');
fprintf('>> bode(sys)\n');
    sys = sss(A,B,C,D,E); clear A B C D E
    disp(sys);
    spysys(sys);
    [~,~,w] = bode(sys);%get frequency range
    figure;bode(sys);%plot
customPause

%   *Modal reduction
fprintf('The first model reduction procedure we would like to test is\n');
fprintf('MODAL REDUCTION, i.e. the truncation of less dominant modes in the\n');
fprintf('system. There are several ways to define "dominance". The default\n');
fprintf('choice in "modalMOR", the MORLab function for modal reduction, is\n');
fprintf('to select the eigenmodes corresponding to the first q eigenvalues \n');
fprintf('with smallest magnitude (q is the desired reduced order). This can \n');
fprintf('be computed quite efficiently using the "eigs" function.\n');
customPause

    q = round(size(sys.A,1)/3);
    isEven = ~(round(q/2)-q/2); if ~isEven, q=q-1; end 
fprintf('For this example, we choose a reduced order of %i, which \n',q);
fprintf('corresponds to approximately a third of the original order.\n');
customPause

fprintf('\nThe reduction is performed by calling:\n');
fprintf('>> sysrModal = modalMOR(sys,q)\n');
    tic, sysrModal = modalMOR(sys,q); tModal = toc;
fprintf('(This reduction took %4.2fs. The results will be shown later.)\n',tModal);
    h2ErrorModal = norm(sys-sysrModal);
    hInfErrorModal = norm(sys-sysrModal,Inf);
customPause
    
%   *Balanced Truncation
fprintf('\nThe next reduction method we wish to test is BALANCED TRUNCATION.\n');
fprintf('This reduction is based on bringing the system to a balanced form,\n');
fprintf('where observability and controllability Gramians are equal and \n');
fprintf('diagonal, and truncating states that are both badly controllable\n');
fprintf('and observable. The reduction order is the same as in the modal \n');
fprintf('case. \n');
customPause

fprintf('The reduction is performed by calling:\n');
fprintf('>> sysrTbr = TBR(sys,q)\n');
tic, sysrTbr = TBR(sys,q); tTbr = toc;
fprintf('(This reduction took %4.2fs. The results will be shown later.)\n',tTbr);
    h2ErrorTbr = norm(sys-sysrTbr);
    hInfErrorTbr = norm(sys-sysrTbr,Inf);
customPause
    
%   *Rational Krylov   
fprintf('\nThe next reduction method we wish to test is implicit moment\n');
fprintf('matching by Krylov subspace methods (RATIONAL KRYLOV). In this \n');
fprintf('procedure, moments, i.e. Taylor series coefficients, of the original\n');
fprintf('model at chosen complex frequencies s0 are matched implicitly by\n');
fprintf('projecting the dynamic equations onto input and output Krylov \n');
fprintf('subspaces.\n');
customPause

fprintf('Since we assume no information about the original model is known \n');
fprintf('we choose to match all moments at the origin. The reduced order q \n');
fprintf('is the same as before, we perform two-sided reduction, hence 2q \n');
fprintf('moments of the original model are matched. \n');
customPause

fprintf('\nThe definition of s0 can be done either by listing all frequenceies \n');
fprintf('in a row vector...\n');
fprintf('>> s0 = [0 0 ... 0]\n');
fprintf('...or in a matrix having frequencies in the 1st and multiplicity in the\n');
fprintf('2nd row.\n');
fprintf('>> s0 = [0;q]\n');
    s0 = [0;q];
customPause

fprintf('The reduction is performed by calling:\n');
fprintf('>> sysrRk = RK(sys,s0,s0)\n');
fprintf('(The first "s0" defines the frequencies for the input, the second "s0"\n');
fprintf('for the output Krylov subspace.)\n');

    tic, sysrRk = RK(sys,s0,s0); tRk = toc;
fprintf('(This reduction took %4.2fs. The results will be shown later.)\n',tRk);
    h2ErrorRk = norm(sys-sysrRk);
    hInfErrorRk = norm(sys-sysrRk,Inf);
customPause
    
%   *IRKA
fprintf('\nAnother algorithm based on rational Krylov which is aimed at \n');
fprintf('finding a set of complex frequencies s0 such that the reduced model\n');
fprintf('is a locally H2 optimal approximation of the original, is the \n');
fprintf('Iterative Rational Krylov Algoritm (IRKA).\n');
customPause

fprintf('To initialize the iteration in IRKA, an initial guess of complex  \n');
fprintf('frequencies s0 has to be provided. Also in this case, we assume \n');
fprintf('no prior knowledge on the dynamics is given and initialize IRKA\n');
fprintf('at the origin:\n');
fprintf('>> s0 = [0 0 ... 0]\n');
customPause

fprintf('The reduction is performed by calling:\n');
fprintf('>> sysrIrka = IRKA(sys,s0)\n');
    tic, sysrIrka = IRKA(sys,s0); tIrka = toc;
fprintf('(This reduction took %4.2fs. The results will be shown later.)\n',tIrka);
    h2ErrorIrka = norm(sys-sysrIrka);
    hInfErrorIrka = norm(sys-sysrIrka,Inf);
customPause
    
%   *Comparison of the results
fprintf('\nCOMPARISON OF THE RESULTS:\n');
fprintf('After reducing the system with several different strategies, we \n');
fprintf('are now ready to compare the reduction results.\n');
fprintf('To this end, we will compare bode plots, error norms and execution\n');
fprintf('times graphically.\n');
customPause

    labs = {'modal','tbr','rk','irka'};

    figure; bode(sys,w,'k-'); hold on
            bode(sysrModal,w,'b--');
            bode(sysrTbr,w,'r-.');
            bode(sysrRk,w,'g--');
            bode(sysrIrka,w,'c-.');
    legend('original',labs{:},'Location','NorthEast');
    
    figure; bar([h2ErrorModal,h2ErrorTbr,h2ErrorRk,h2ErrorIrka]);
            title('H2 error');
            set(gca,'XTIckLabel',labs);
            
    figure; bar([hInfErrorModal,hInfErrorTbr,hInfErrorRk,hInfErrorIrka]);
            title('HInf error');
            set(gca,'XTIckLabel',labs);
figure; bar([tModal,tTbr,tRk,tIrka]);
            title('Execution times /s');
            set(gca,'XTIckLabel',labs);
customPause


rule
demoEnd
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
            isValid = 0;
        otherwise    
            warning('Your input was not valid')
            fprintf(['\nChoose one of the following benchmark models by pressing...\n',...
                '\t "1" \t ..for the Los Angeles Hospital model (N=48)\n',...
                '\t "2" \t ..for a 1D-beam model (N=348)\n']);
            modelChosen = input('Model choice: ');
    end
end