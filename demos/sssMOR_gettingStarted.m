function sssMOR_gettingStarted
% SSSMOR_GETTINGSTARTED - Introductory demo to sssMOR toolbox 
% 
% Syntax:
%       SSSMOR_GETTINGSTARTED
% 
% Description:
%       This function can be executed as it is and will guide the 
%       user in the command window through an introductory journey in the 
%       capabilities and advantages of using the model reduction functions.
%
% See Also:
%       sss_gettingStarted
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State-Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
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
fprintf(['\t sssMOR - A sparse state-space and model reduction toolbox\n']); 
rule
fprintf(['\tDeveloped at the Chair of Automatic Control, TUM\n']);
fprintf(['\tNote: for academic use only.\n']);
rule
%%  Model Order Reduction
clear

fprintf(['\t Model Order Reduction (MOR)\n\n']);

fprintf('Besides allowing the definition of sparse dynamic systems, the \n');
fprintf('main goal of the sssMOR toolbox is to provide classic and state-of- \n');
fprintf('the-art MOR functions for academic research and teaching.\n');

fprintf('The functions are implemented to work in combination with the sss- \n');
fprintf('class and comprise classical MOR techniques, such as modal reduction, \n');
fprintf('balanced truncation or rational Krylov, as well as more advanced \n');
fprintf('techniques developed in recent years. \n');
customPause

fprintf('\nWe will now demonstrate the main MOR capabilities of the sssMOR \n');
fprintf('toolbox through an example.\n');

%   *Selection of a benchmark model
[sysName,A,B,C,D,E] = selectModel;

%   *Display some information on the model
fprintf('\nYou chose the "%s" model. \n',sysName);

fprintf('\nDefine a sss-object and show some information about it.\n');
fprintf('>> sys = sss(A,B,C)\n');
fprintf('>> disp(sys)\n');
fprintf('>> spy(sys)\n');
fprintf('>> bode(sys)\n');
    sys = sss(A,B,C); clear A B C D E
    disp(sys);
    spy(sys);
    [~,~,w] = bode(sys);%get frequency range
    figure;bode(sys);%plot
customPause

fprintf(['Only the first output and the first input of the CDplayer will\n'...
        'be taken into consideration in the reduction.\n']);
fprintf(['This can be done using the function truncate:\n'...
        '>>sys=truncate(sys,1,1);\n'])
    sys=truncate(sys,1,1);
    fprintf('>>disp(sys)\n');
    disp(sys);
customPause
%   *Modal reduction
fprintf('The first model reduction procedure we would like to test is\n');
fprintf('MODAL REDUCTION, i.e. the truncation of less dominant modes in the\n');
fprintf('system. There are several ways to define "dominance". The default\n');
fprintf('choice in "modalMor", the sssMOR function for modal reduction, is\n');
fprintf('to select the eigenmodes corresponding to the first q eigenvalues \n');
fprintf('with smallest magnitude (q is the desired reduced order). This can \n');
fprintf('be computed quite efficiently using the "eigs" function.\n');
customPause
    if strcmp(sysName,'CDplayer')
    q = round(size(sys.A,1)/6);
    isEven = ~(round(q/2)-q/2); if ~isEven, q=q-1; end 
fprintf('For this example, we choose a reduced order of %i, which \n',q);
fprintf('corresponds to approximately a one-sixth of the original order.\n');
    else
    q = round(size(sys.A,1)/3);
    isEven = ~(round(q/2)-q/2); if ~isEven, q=q-1; end 
fprintf('For this example, we choose a reduced order of %i, which \n',q);
fprintf('corresponds to approximately a third of the original order.\n');
    end
customPause

fprintf('\nThe reduction is performed by calling:\n');
fprintf('>> sysrModal = modalMor(sys,q)\n');
    Opts.real = 'real';
    %modalMor can fail with building, q=16 to find matching eigenvalues at the
    %first try (eigs in initialized randomly)
    try
        tic, sysrModal = modalMor(sys,q,Opts); tModal = toc;
    catch
        tic, sysrModal = modalMor(sys,q,Opts); tModal = toc;
    end
fprintf('(This reduction took %4.2fs. The results will be shown later.)\n',tModal);
    h2ErrorModal = norm(sys-sysrModal);
    hInfErrorModal = norm(sys-sysrModal,Inf);
customPause
    
%   *Balanced Truncation
fprintf('\nThe next reduction method we wish to test is BALANCED TRUNCATION.\n');
fprintf('This reduction is based on bringing the system to a balanced form,\n');
fprintf('where observability and controllability Gramians are equal and \n');
fprintf('diagonal, and truncating states variables that are both badly contro-\n');
fprintf('llable and observable. The reduction order is the same as in the modal \n');
fprintf('case. \n');
customPause

fprintf('The reduction is performed by calling:\n');
fprintf('>> sysrTbr = tbr(sys,q)\n');
tic, sysrTbr = tbr(sys,q); tTbr = toc;
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

fprintf('\nThe definition of s0 can be done either by listing all frequencies \n');
fprintf('in a row vector...\n');
fprintf('>> s0 = [0 0 ... 0]\n');
fprintf('...or in a matrix having frequencies in the 1st and multiplicity in the\n');
fprintf('2nd row.\n');
fprintf('>> s0 = [0;q]\n');
    s0 = [0;q];
customPause

fprintf('The reduction is performed by calling:\n');
fprintf('>> sysrRk = rk(sys,s0,s0)\n');
fprintf('(The first "s0" defines the frequencies for the input, the second "s0"\n');
fprintf('for the output Krylov subspace.)\n');

    tic, sysrRk = rk(sys,s0,s0); tRk = toc;
fprintf('(This reduction took %4.2fs. The results will be shown later.)\n',tRk);
    h2ErrorRk = norm(sys-sysrRk);
    hInfErrorRk = norm(sys-sysrRk,Inf);
customPause
    
%   *IRKA
fprintf('\nAnother algorithm based on rational Krylov which is aimed at \n');
fprintf('finding a set of complex frequencies s0 such that the reduced model\n');
fprintf('is a locally H2-optimal approximation of the original, is the \n');
fprintf('Iterative Rational Krylov Algoritm (IRKA).\n');
customPause

fprintf('To initialize the iteration in IRKA, an initial guess of complex  \n');
fprintf('frequencies s0 has to be provided. Also in this case, we assume \n');
fprintf('no prior knowledge on the dynamics is given and initialize IRKA\n');
fprintf('at the origin:\n');
fprintf('>> s0 = [0 0 ... 0]\n');
customPause

fprintf('The reduction is performed by calling:\n');
fprintf('>> sysrIrka = irka(sys,s0)\n');
    tic, sysrIrka = irka(sys,s0); tIrka = toc;
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
demoEndsssMOR
close all
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
                '\t "2" \t ..for a CD-Player model (N=120)\n']);
modelChosen = input('Model choice: ');

if isempty(modelChosen),modelChosen = 0;end
isValid = 0;
while ~isValid
    switch modelChosen
        case 1 %building
            sysName = 'building';
            LoadData = load('building.mat'); %assumes it is in path
            
            A = LoadData.A; B = LoadData.B; C = LoadData.C;
            E = speye(size(A));
            D = 0;
            isValid = 1;
        case 2 %CDplayer
            sysName = 'CDplayer';
            LoadData = load('CDplayer.mat');%assumes it is in path
            
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
                '\t "2" \t ..for a CD-Player model (N=120)\n']);
            modelChosen = input('Model choice: ');
            if isempty(modelChosen),modelChosen = 0;end
    end
end

